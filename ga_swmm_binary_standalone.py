import numpy as np
import random
from typing import Tuple, List
import matplotlib.pyplot as plt
import re
#import chardet
import pandas as pd
from swmmtoolbox import swmmtoolbox
from pyswmm import Simulation, Nodes
import copy
import time
import networkx as nx
import subprocess

# =============== User Parameters (you can edit) ===============
# Economics
#om_lid_bc = 0.08        # bioretention O&M rate (per year relative to capex)
#om_lid_pp = 0.04        # permeable pavement O&M rate
#discount_rate = 0.02    # discount rate
#project_lifetime = 30   # years


# GA hyperparameters
population_size = 50   # 种群规模 (use 3/5/8 to test)
#num_genes = 329        # 每个个体的基因长度
num_generations = 50    # 最大迭代次数
crossover_rate = 0.8    # probability of crossover
mutation_rate  = 0.006   # per-bit flip probability
#ELITE_K        = 2
#TOURN_K        = 3

# Problem dimensions fixed by your mapping
NUM_PIPES        = 83
#NUM_OUTLET_CANDS = 8
NUM_SUBCATCH     = 18

# Normalization bounds
LCC_MIN, LCC_MAX = 3e4, 3e6
Q_MAX_MIN, Q_MAX_MAX = 0, 3e4
EPS = 1e-6

# Objective weights (minimize w1*LCC_norm + w2*Qmax_norm)
w1 = 0.5  # weight for LCC
w2 = 0.5  # weight for max outflow

# 2-bit diameter code -> commercial diameters (meters). Customize if needed.
# Load conduits
conduits_info = np.loadtxt("[CONDUITS].txt")
# Load commercial diameters
commercial_df = pd.read_csv("Commercial.txt", sep=r"\s+", header=None)
commercial_df.columns = ["Diameter_m", "UnitCost"]
commercial_diameters = commercial_df["Diameter_m"].tolist()
# cost dictionary from file (if needed)
cost_per_diameter = dict(zip(commercial_df["Diameter_m"], commercial_df["UnitCost"]))
# 离散可选权重（管径）列表，严格递增
weight_list = commercial_diameters.copy()
weight_to_state = { w: weight_list.index(w) + 1 for w in weight_list } # 1-based index

# Pipe lengths (meters): simple synthetic dataset, per-pipe. Replace with your real data if available.
#_rng = np.random.default_rng(42)
#PIPE_LENGTHS_M = _rng.uniform(50.0, 150.0, size=NUM_PIPES).astype(float)
loops = np.loadtxt("Loops.txt", dtype=int)
# Subcatchment areas (m^2), 18 values; replace with real data if available.
subcatch_info    = np.loadtxt("[SUBCATCHMENTS].txt", dtype=float)
catchment_areas_ha = subcatch_info[:, 3] # to calculate the LID areas
catchment_areas = catchment_areas_ha * 10000  # km² → m²
subarea_info = np.loadtxt("[SUBAREAS].txt", dtype=float)
# read Subcatchout.txt
# each row：[subcatchment_ID, OutletNode]
subcatchout_info = np.loadtxt("Subcatchout.txt", dtype=int, delimiter=",")


# pipe_ids 也要手工指定
pipe_ids = [2,17,12,10,8,15,16,23,24,30,32,37,44,55,41,48,53,50]
n_sc = subcatch_info.shape[0]
# create a dict: {ID → OutletNode}
#outlet_map = { row[0]: row[1] for row in subcatchout_info }
coords_info = np.loadtxt("[COORDINATES].txt", dtype=float, delimiter=",").reshape(-1, 3)

junctions_info = np.loadtxt("[JUNCTIONS].txt", dtype=float)
elev_map = { int(row[0]): float(row[1]) for row in junctions_info }

all_nodes = np.unique(conduits_info[:,1:3].astype(int)) # number of nodes: 48
# node_list 保证顺序稳定
node_list = sorted(all_nodes)
node_to_idx = { node:i for i,node in enumerate(node_list) }
M = len(node_list)  # 出口槽数（例如 48）

def static_split(conduits_info, coords_info, outlet_nodes, pipe_ids):
    """
    在 conduits_info 上静态地把 outlet_nodes[i] 插入到 pipe pipe_ids[i] 上，
    产生一个“master” conduits_info 与对应 coords_info。
    要求 outlet_nodes, pipe_ids 长度相同。
    返回 master_coords_info, master_conduits_info。
    """
    coords_dict = {int(r[0]): (r[1], r[2]) for r in coords_info}
    master_nodes = coords_info.copy().tolist()
    master_conduits = [row.copy() for row in conduits_info]

    next_node_id = int(coords_info[:,0].max()) + 1
    next_pipe_id = int(conduits_info[:,0].max()) + 1

    for outlet_node, orig_pipe in zip(outlet_nodes, pipe_ids):
        # （1）在 master_conduits 里找出 orig_pipe 的那条管道及其索引
        for idx, row in enumerate(master_conduits):
            if int(row[0]) == orig_pipe:
                break
        else:
            raise RuntimeError(f"找不到管道 ID={orig_pipe}")
        u, v = int(row[1]), int(row[2])
        xu, yu = coords_dict[u]; xv, yv = coords_dict[v]
        xo, yo = coords_dict[outlet_node]

        # 2) 插值新节点
        # 如果 outlet_node 已经在 coords_dict 中，就复用它，不新建
        if outlet_node in coords_dict:
            new_id = outlet_node
        else:
            new_id = next_node_id
            next_node_id += 1
            master_nodes.append([new_id, xo, yo])

        # 3) 把管道 u->v 拆成 u->new_id 和 new_id->v
        #    用原 pipe id 复用第一段，用新 pipe id 做第二段
        row1 = row.copy()
        row1[2] = new_id  # u->new_id
        row2 = row.copy()
        row2[0] = next_pipe_id; next_pipe_id += 1
        row2[1] = new_id  # new_id->v

        # （4）用 pop(idx) 取出原来的 row，再在相同位置插入两段
        master_conduits.pop(idx)
        master_conduits.insert(idx, row1)
        master_conduits.insert(idx+1, row2)

    return np.array(master_nodes), np.array(master_conduits)

# 1) 一次性做静态切割，得到 master network
master_coords, master_conduits = static_split(
    conduits_info,
    coords_info,
    outlet_nodes=subcatchout_info[:,1].tolist(),
    pipe_ids=pipe_ids
)

# 2) 覆盖原有
coords_info = master_coords
conduits_info = master_conduits

# 更新全局节点字典
coords_dict = {int(r[0]):(r[1],r[2]) for r in coords_info}

# Rebuild your edge→index map (in case it changed)
edge_to_pipe_idx = {
    (int(row[1]), int(row[2])): i
    for i, row in enumerate(conduits_info)
}
N_pipes = conduits_info.shape[0] # 83
# G_full is your directed “master graph”
G_full = nx.DiGraph()
for (u, v), idx in edge_to_pipe_idx.items():
    G_full.add_edge(u, v)

# build an undirected copy with a `weight` attribute on each edge
G_und = G_full.to_undirected()
for (u, v), idx in edge_to_pipe_idx.items():
    G_und[u][v]['idx']    = idx
    G_und[u][v]['weight'] = float(conduits_info[idx, 3])  # use length as weight

# Compute the MST over that graph
T = nx.minimum_spanning_tree(G_und, weight='weight') # undirected 
# print T edge list without weights

print("MST edges (undirected):", T.edges(data=False))
# size of the MST
print("MST size (edges):", T.number_of_edges())


# Extract the indices of the edges in the MST
mst_indices = { data['idx'] for u, v, data in T.edges(data=True) }

# Now you have the correct MST for the 83-link network
print("Recomputed MST covers pipe indices:", sorted(mst_indices))

# ─── NOW rebuild the skeleton DAG on the new conduits_info ───
pipe_skel_dag = nx.DiGraph()


# 1) add every MST pipe‐index as a node (index)
for idx in mst_indices:
    pipe_skel_dag.add_node(idx)
# 2) for each MST edge idx→jdx if idx.to_node == jdx.from_node (edge)
for idx in mst_indices:
    to_node = int(conduits_info[idx, 2])
    for jdx in mst_indices:
        if int(conduits_info[jdx, 1]) == to_node:
            pipe_skel_dag.add_edge(idx, jdx)
# print the skeleton DAG with edges
print("Skeleton DAG edges:", pipe_skel_dag.edges())


# 3) topological order and predecessors map
mst_pipe_order = list(nx.topological_sort(pipe_skel_dag))
skel_preds     = { p: list(pipe_skel_dag.predecessors(p))
                   for p in pipe_skel_dag.nodes() }

# 子流域出口集合
subcatch_outlets = set(subcatchout_info[:,1].tolist())

# 找所有 G_full 中出度为 0 的节点
candidates = [n for n in G_full.nodes() if G_full.out_degree(n)==0]

# 排除子流域接口节点，只保留真正管网末端
system_outfalls = [n for n in candidates if n not in subcatch_outlets]
if not system_outfalls:
    raise RuntimeError("No valid system outfall found!")
# 随机选两个作为出水口
system_outfall_node = random.sample(system_outfalls, 2)
OUTLET_FIXED_MASK = np.zeros(M, dtype=np.uint8)
for n in system_outfall_node:
    OUTLET_FIXED_MASK[node_to_idx[n]] = 1

# =============== Genome utilities (exact 329-bit mapping) ===============
# Index ranges (inclusive) for readability
IDX_PIPE_EN_START = 0
IDX_PIPE_EN_END   = IDX_PIPE_EN_START + NUM_PIPES - 1                  # 0..82

IDX_OUTLET_START  = IDX_PIPE_EN_END + 1                                # 83
IDX_OUTLET_END    = IDX_OUTLET_START + M - 1                           # 83..(82+M)

IDX_DIAM_START    = IDX_OUTLET_END + 1                                 # ...
IDX_DIAM_END      = IDX_DIAM_START + (NUM_PIPES*4) - 1                    

IDX_BC_START      = IDX_DIAM_END + 1                                   
IDX_BC_END        = IDX_BC_START + (NUM_SUBCATCH*2) - 1                 

IDX_PP_START      = IDX_BC_END + 1                                     
IDX_PP_END        = IDX_PP_START + (NUM_SUBCATCH*2) - 1                 

TOTAL_BITS = IDX_PP_END + 1    
# 打印各段区间
print(f"Genome segments (0-indexed, inclusive):")
print(f"PIPE_EN: {IDX_PIPE_EN_START}-{IDX_PIPE_EN_END}")
print(f"OUTLET_M: {IDX_OUTLET_START}-{IDX_OUTLET_END}")
print(f"DIAM_4B: {IDX_DIAM_START}-{IDX_DIAM_END}")
print(f"BC_BITS: {IDX_BC_START}-{IDX_BC_END}")
print(f"PP_BITS: {IDX_PP_START}-{IDX_PP_END}")

#assert TOTAL_BITS == 321+M, f"Genome must be 329 bits, got {TOTAL_BITS}"
# ===== Sanity: 片段长度与不重叠 =====
assert IDX_PIPE_EN_END - IDX_PIPE_EN_START + 1 == NUM_PIPES
assert IDX_OUTLET_END  - IDX_OUTLET_START  + 1 == M
assert IDX_DIAM_END    - IDX_DIAM_START    + 1 == NUM_PIPES * 4  # 4位/管
assert IDX_BC_END      - IDX_BC_START      + 1 == 2 * NUM_SUBCATCH
assert IDX_PP_END      - IDX_PP_START      + 1 == 2 * NUM_SUBCATCH
assert TOTAL_BITS == IDX_PP_END + 1

# 不重叠（含式区间）
segs = [
    (IDX_PIPE_EN_START, IDX_PIPE_EN_END, "PIPE_EN"),
    (IDX_OUTLET_START,  IDX_OUTLET_END,  "OUTLET_M"),
    (IDX_DIAM_START,    IDX_DIAM_END,    "DIAM_4B"),
    (IDX_BC_START,      IDX_BC_END,      "BC_BITS"),
    (IDX_PP_START,      IDX_PP_END,      "PP_BITS"),
]
segs_sorted = sorted(segs, key=lambda t: t[0])
for (a1,b1,n1),(a2,b2,n2) in zip(segs_sorted, segs_sorted[1:]):
    assert b1 < a2, f"Overlap: {n1}({a1}-{b1}) vs {n2}({a2}-{b2})"
print(f"[GENOME OK] total={TOTAL_BITS} :: " + " | ".join(f"{n}:{a}-{b}" for a,b,n in segs_sorted))

def gray_to_bin(g):
    x = g
    while g:
        g >>= 1
        x ^= g
    return x

def decode_pipe_diams_4bit(g_bits, commercial_diameters):
    #diams_idx = []
    diams = np.zeros(NUM_PIPES, dtype=float)
    for i in range(NUM_PIPES):
        off = IDX_DIAM_START + 4*i
        # 小端读取（保持你原有风格 b0 最低位）
        code = (int(g_bits[off+0]) |
               (int(g_bits[off+1])<<1) |
               (int(g_bits[off+2])<<2) |
               (int(g_bits[off+3])<<3))
        idx = gray_to_bin(code)            # 0..15
        idx = min(idx, len(commercial_diameters)-1, 9)  # 夹到 0..9
        diams[i] = float(commercial_diameters[idx])
        #diams_idx.append(idx)
    #print("diam idx hist:", np.bincount(diams_idx, minlength=10))
    return diams


def decode_genome_M(g: np.ndarray, node_list:list):
    """→ pipe_enabled(83,), pipe_diams_m(83,), outlet_nodes(list[0..7]), lid_vars(36,)"""
    assert isinstance(g, np.ndarray) and g.ndim == 1 and g.size == TOTAL_BITS, \
        f"genome length {getattr(g,'size',None)} != TOTAL_BITS={TOTAL_BITS}"
    # 1) 管道启用
    pipe_enabled = g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1].astype(np.int32)

    # 2) 出口（M 位）：把置 1 的位映射回 node_list 的节点ID
    outlet_bits = g[IDX_OUTLET_START:IDX_OUTLET_END+1]
    outlet_nodes = [node_list[i] for i, b in enumerate(outlet_bits) if b > 0]

    # 3) 管径（2位/管）
    '''
    d_bits = g[IDX_DIAM_START:IDX_DIAM_END+1]
    pipe_diameters = np.zeros(NUM_PIPES, dtype=float)
    for i in range(NUM_PIPES):
        b0 = int(d_bits[2*i + 0]); b1 = int(d_bits[2*i + 1])
        code = (b1 << 1) | b0    # 小端
        pipe_diameters[i] = float(commercial_diameters[code])
    '''
    pipe_diameters = decode_pipe_diams_4bit(g, commercial_diameters)  
    # 4) LID：BC/PP 各 2 位/子汇：比例: 00,01,10,11 → 0, 0.025, 0.05, 0.10
    lid_vars = np.zeros(2*NUM_SUBCATCH, dtype=float)  # 偶数位BC，奇数位PP（与你的check_lid_constraints对齐）
    def _decode_lid(start_idx: int, offset: int):
        bits = g[start_idx:start_idx + 2*NUM_SUBCATCH]
        for i in range(NUM_SUBCATCH):
            _LID_LEVELS = (0.0, 0.025, 0.05, 0.10)
            code = (int(bits[2*i]) << 1) | int(bits[2*i + 1])   # 00/01/10/11 → 0/1/2/3
            lid_vars[2*i + offset] = _LID_LEVELS[code]
    _decode_lid(IDX_BC_START, 0)  # 偶数位
    _decode_lid(IDX_PP_START, 1)  # 奇数位
    return pipe_enabled, pipe_diameters, outlet_nodes, lid_vars

def check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes,
                             edge_to_pipe_idx, G_full,
                             subcatch_outlets=None,
                             require_single_component=True,     # 整网一个弱连通组件
                             require_reachability=False          # 每个节点都能到至少一个出口
                             ):

    if not outlet_nodes:                  # 至少有出口
        return False

    # 启用子图（有向）
    enabled_edges = [(u,v) for (u,v),idx in edge_to_pipe_idx.items() if pipe_enabled[idx]]
    if not enabled_edges:
        return False
    H = G_full.edge_subgraph(enabled_edges).copy()
    U = H.to_undirected()

    # 出口必须在子图中，并且在 H 里是末端（出度=0），且不属于子流域接口
    for o in outlet_nodes:
        if o not in H: return False
        if H.out_degree(o) != 0: return False
        if subcatch_outlets and o in subcatch_outlets: return False

    # 弱连通性：整网一个组件，或放宽为“每组件至少有一个出口”
    if require_single_component:
        if not nx.is_connected(U):
            return False
    else:
        for comp_nodes in nx.connected_components(U):
            if not any(o in comp_nodes for o in outlet_nodes):
                return False

    # 上下游管径单调：上游 ≤ 下游（仅对启用边）
    for (u,v), idx_d in edge_to_pipe_idx.items():
        if not pipe_enabled[idx_d]:
            continue
        d_down = pipe_diameters[idx_d]
        for (w,uu), idx_u in edge_to_pipe_idx.items():
            if uu == u and pipe_enabled[idx_u]:
                d_up = pipe_diameters[idx_u]
                if d_up > 0 and d_down > 0 and d_up > d_down + 1e-9:
                    return False

    # 可达性检查：从所有出口出发在“反向图”里做多源搜索
    if require_reachability:
        Hr = H.reverse(copy=False)
        seen = set(outlet_nodes)
        for s in outlet_nodes:
            seen |= nx.descendants(Hr, s) | {s}   # 能到达 s 的点（在原图里能流向 s）
        # 只对出现在启用边中的节点要求可达
        if not set(H.nodes()).issubset(seen):
            return False

    return True

## LID constraints
def check_lid_constraints(lid_vars, catchment_areas, max_frac=0.1):
    """
    lid_vars: 1-D array of length 36,
       even indices are BC fractions, odd indices are PP fractions.
    catchment_areas: array of length 18 (areas in m²).
    max_frac: maximum total fraction (e.g. 0.1 for 10%).
    """
    n_sc = len(catchment_areas)    # should be 18
    # split lid_vars into pairs
    bc_fracs = lid_vars[0::2]      # BC for each catchment
    pp_fracs = lid_vars[1::2]      # PP for each catchment
    for i in range(n_sc):
        # convert fraction back to area
        total_frac = bc_fracs[i] + pp_fracs[i]
        if total_frac > max_frac + 1e-6:  # allow tiny epsilon
            return False
    return True

def init_population_feasible(population_size: int,
                             add_prob: float = 0.15, # or 0.5
                             tries_per_ind: int = 30,
                             seed: int = 2025) -> np.ndarray:
    """
    产出形状 (population_size, TOTAL_BITS) 的初始种群（M 位出口）。
    构造策略：骨架连通 + 随机加边、随机管径后单调化、从候选末端挑2个出口、LID在可行域内采样。
    生成后用你已有的 check_* 再验一遍；失败重试，最后兜底。
    依赖全局：G_full, edge_to_pipe_idx, mst_indices, node_list,
             subcatch_outlets, system_outfall_node, catchment_areas。
    """
    # 用 2bit 表示 4 档：00/01/10/11 → 0/2.5%/5%/10%
    _LID_LEVELS = [0.0, 0.025, 0.05, 0.10]
    def _encode_four_bits(g, off, idx, use_gray=True):
        code = bin_to_gray(idx) if use_gray else idx
        for k in range(4):
            g[off+k] = (code >> k) & 1

    def _encode_lid_pair(dst: np.ndarray, off: int, frac: float):
        if frac < 0.0125: code = 0
        elif frac < 0.0375: code = 1
        elif frac < 0.075:  code = 2
        else:               code = 3
        en = code >> 1      # 高位
        ab = code & 1       # 低位
        dst[off + 0] = en
        dst[off + 1] = ab
    np.random.seed(seed); random.seed(seed)
    pop = np.zeros((population_size, TOTAL_BITS), dtype=np.uint8)
    # MST 骨架 → 保证连通
    backbone = np.zeros(NUM_PIPES, dtype=np.uint8)
    for idx in mst_indices:
        backbone[idx] = 1

    for k in range(population_size):
        for _try in range(tries_per_ind):
            g = np.zeros(TOTAL_BITS, dtype=np.uint8)

            # 1) 管道启用：骨架 + 随机加边
            enabled = backbone.copy()
            flip = (np.random.rand(NUM_PIPES) < add_prob)
            enabled = np.where((enabled == 1) | (flip == 1), 1, 0).astype(np.uint8)
            g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1] = enabled

            # 2) 管径 2bit×83：随机档位(0..3) + 单调化
            codes = np.random.randint(0, min(10, len(commercial_diameters)), size=NUM_PIPES, dtype=int)
            for _ in range(6):  # limited iteration to avoid too monotonic
                changed = False
                for (u, v), idx_down in edge_to_pipe_idx.items():
                    if enabled[idx_down] == 0: 
                        continue
                    for (w, uu), idx_up in edge_to_pipe_idx.items():
                        if uu == u and enabled[idx_up] and codes[idx_up] > codes[idx_down]:
                            codes[idx_down] = codes[idx_up]
                            changed = True
                if not changed: break
            # 写回4位（小端）
            for i in range(NUM_PIPES):
                _encode_four_bits(g, IDX_DIAM_START + 4*i, int(codes[i]), use_gray=True)
            # 3) 出口 M 位：从“启用子图”的末端（出度=0，排除子流域接口）挑 2 个
            '''
            H = nx.DiGraph()
            H.add_nodes_from(G_full.nodes())
            H.add_edges_from([(u, v) for (u, v), idx in edge_to_pipe_idx.items() if enabled[idx]])
            sinks = [n for n in H.nodes() if H.out_degree(n) == 0 and n not in subcatch_outlets]

            # 选择两个（优先末端，不足从 node_list 补齐）
            #chosen = [n for n in system_outfall_node if n in sinks]  # 与启用子图一致的
            chosen = []
            if len(sinks) >= 2:
                # 你也可以改成 random.sample(sinks, 2)
                chosen = sinks[:2]
            else:
                chosen = sinks.copy()
                for n in node_list:
                    if n not in chosen:
                        chosen.append(n)
                    if len(chosen) == 2: break

            outM = np.zeros(M, dtype=np.uint8)
            for n in chosen[:2]:
                outM[node_to_idx[n]] = 1
            '''
            outM = OUTLET_FIXED_MASK.copy()
            g[IDX_OUTLET_START:IDX_OUTLET_END+1] = outM  # —— 关键：写入 M 位出口段

            # 4) LID：每子汇 4 档采样（00/01/10/11），并修复 bc+pp ≤ 10%
            for i in range(NUM_SUBCATCH):
                # 随机选 code ∈ {0,1,2,3}
                code_bc = np.random.randint(0, 4)
                code_pp = np.random.randint(0, 4)
                bc = _LID_LEVELS[code_bc]
                pp = _LID_LEVELS[code_pp]

                # 若超 10%，保留较大者，较小者下调到 ≤ (0.10 - 大者) 的最大可行档
                if bc + pp > 0.10 + 1e-12:
                    if bc >= pp:
                        # 给 PP 选一个不超过剩余额度的最大档
                        remain = max(0.0, 0.10 - bc)
                        pp = max([v for v in _LID_LEVELS if v <= remain] or [0.0])
                    else:
                        remain = max(0.0, 0.10 - pp)
                        bc = max([v for v in _LID_LEVELS if v <= remain] or [0.0])

                _encode_lid_pair(g, IDX_BC_START + 2*i, bc)
                _encode_lid_pair(g, IDX_PP_START + 2*i, pp)

            # 5) 用你已有的解码+约束函数验一遍（保持与你现有实现一致）
            pipe_enabled, pipe_diameters, outlet_nodes, lid_vars = decode_genome_M(g, node_list)
            try:
                ok_struct = check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes,edge_to_pipe_idx=edge_to_pipe_idx, G_full=G_full, subcatch_outlets=subcatch_outlets)
            except Exception:
                ok_struct = False
            ok_lid = check_lid_constraints(lid_vars, catchment_areas, max_frac=0.10)

            if ok_struct and ok_lid:
                pop[k] = g
                break
        else:
            # 兜底：只开骨架、最小管径、前两个出口、LID=0
            g = np.zeros(TOTAL_BITS, dtype=np.uint8)
            g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1] = backbone
            for i in range(NUM_PIPES):
                off = IDX_DIAM_START + 4*i
                g[off + 0] = 0; g[off + 1] = 0; g[off + 2] = 0; g[off + 3] = 0
            outM = np.zeros(M, dtype=np.uint8); outM[0:2] = 1
            g[IDX_OUTLET_START:IDX_OUTLET_END+1] = outM
            pop[k] = g
    return pop

# ====== 1) LID 约束解释 ======
def explain_lid(lid_vars, catchment_areas, max_frac=0.10, max_show=10):
    bc = lid_vars[0::2]; pp = lid_vars[1::2]
    bad = []
    for i, (b, p) in enumerate(zip(bc, pp)):
        s = b + p
        if s > max_frac + 1e-6:
            bad.append((i, b, p, s))
    if bad:
        print(f"[LID] Violations: {len(bad)} subcatchments exceed {max_frac:.2f}. Showing up to {max_show}:")
        for i, b, p, s in bad[:max_show]:
            print(f"   - Subcatch #{i}: BC={b:.3f}, PP={p:.3f}, SUM={s:.3f}")
        return False
    print("[LID] OK: all subcatchments satisfy BC+PP <= 0.10")
    return True

# ====== 2) 结构约束解释（带原因）======
def explain_structure(pipe_enabled, pipe_diameters, outlet_nodes,
                      edge_to_pipe_idx, G_full,
                      subcatch_outlets=None,
                      require_single_component=True,
                      require_reachability=False,
                      max_show=10):
    import networkx as nx
    # 0) 出口
    if not outlet_nodes:
        print("[STRUCT] FAIL: no outlet nodes selected.")
        return False
    print(f"[STRUCT] Outlets: {outlet_nodes}")

    # 1) 启用子图
    enabled_edges = [(u,v) for (u,v), idx in edge_to_pipe_idx.items() if pipe_enabled[idx]]
    if not enabled_edges:
        print("[STRUCT] FAIL: no enabled edges.")
        return False
    H = G_full.edge_subgraph(enabled_edges).copy()
    U = H.to_undirected()
    print(f"[STRUCT] Enabled edges: {len(enabled_edges)}, nodes in enabled subgraph: {H.number_of_nodes()}")

    # 2) 出口必须在子图里&是末端&不是子流域接口
    fail = False
    for o in outlet_nodes:
        if o not in H:
            print(f"[STRUCT] FAIL: outlet node {o} NOT in enabled subgraph.")
            fail = True
            continue
        if H.out_degree(o) != 0:
            print(f"[STRUCT] FAIL: outlet node {o} has out_degree={H.out_degree(o)} (not a sink).")
            fail = True
        if subcatch_outlets and (o in subcatch_outlets):
            print(f"[STRUCT] FAIL: outlet node {o} is a subcatch interface (not allowed).")
            fail = True
    if fail: return False

    # 3) 连通性
    if U.number_of_nodes() == 0 or not nx.is_connected(U):
        comps = list(nx.connected_components(U))
        print(f"[STRUCT] FAIL: enabled subgraph has {len(comps)} components (require single component={require_single_component}).")
        if not require_single_component:
            # 放宽规则：每个组件至少一个出口
            ok_all = True
            for k, comp in enumerate(comps):
                has_out = any(o in comp for o in outlet_nodes)
                print(f"   - component {k} size={len(comp)} has_outlet={has_out}")
                if not has_out: ok_all = False
            if ok_all:
                print("[STRUCT] PASS (relaxed): every component has >=1 outlet.")
            else:
                return False
        else:
            return False
    else:
        print("[STRUCT] OK: enabled subgraph is connected.")

    # 4) 上下游直径单调
    mono_bad = []
    for (u,v), idx_d in edge_to_pipe_idx.items():
        if not pipe_enabled[idx_d]: continue
        d_down = pipe_diameters[idx_d]
        for (w,uu), idx_u in edge_to_pipe_idx.items():
            if uu == u and pipe_enabled[idx_u]:
                d_up = pipe_diameters[idx_u]
                if d_up > 0 and d_down > 0 and d_up > d_down + 1e-9:
                    mono_bad.append(((w,uu),(u,v), d_up, d_down))
    if mono_bad:
        print(f"[STRUCT] FAIL: diameter monotonicity violated on {len(mono_bad)} places. Showing up to {max_show}:")
        for ((w,uu),(u,v), up, dn) in mono_bad[:max_show]:
            print(f"   - {w}->{uu} (up={up:.3f}) > {u}->{v} (down={dn:.3f})")
        return False
    print("[STRUCT] OK: upstream diameter <= downstream diameter.")

    # 5) 可达性：每个启用子图节点能流到任一出口
    if require_reachability:
        Hr = H.reverse(copy=False)
        reachable = set()
        for s in outlet_nodes:
            reachable |= nx.descendants(Hr, s) | {s}   # 能到达 s 的点（原图方向能流向 s）
        not_cov = set(H.nodes()) - reachable
        if not_cov:
            print(f"[STRUCT] FAIL: {len(not_cov)} nodes cannot reach any outlet. Showing up to {max_show}:")
            for n in list(not_cov)[:max_show]:
                print(f"   - node {n}")
            return False
        print("[STRUCT] OK: every node can reach an outlet.")
    return True

# ====== 3) 综合调试入口 ======
def debug_constraints_for_individual(individual,
                                     node_list,
                                     catchment_areas,
                                     edge_to_pipe_idx, G_full,
                                     subcatch_outlets=None,
                                     require_single_component=True,
                                     require_reachability=False):
    print("="*60)
    print("Debug constraints for this individual")
    pipe_enabled, pipe_diameters, outlet_nodes, lid_vars = decode_genome_M(individual, node_list)
    print(f"[DECODE] enabled={int(np.sum(pipe_enabled))}/{len(pipe_enabled)} pipes | outlets={outlet_nodes} | "
          f"LID-mean BC={lid_vars[0::2].mean():.3f}, PP={lid_vars[1::2].mean():.3f}")

    ok_lid = explain_lid(lid_vars, catchment_areas, max_frac=0.10)
    ok_struct = explain_structure(pipe_enabled, pipe_diameters, outlet_nodes,
                                  edge_to_pipe_idx=edge_to_pipe_idx, G_full=G_full,
                                  subcatch_outlets=subcatch_outlets,
                                  require_single_component=require_single_component,
                                  require_reachability=require_reachability)
    print(f"[RESULT] LID={ok_lid} | STRUCT={ok_struct}")
    print("="*60)
    return ok_lid and ok_struct

def _bits_to_frac(en, ab):
    return 0.0 if en == 0 else (0.05 if ab == 0 else 0.10)

def _frac_to_bits(frac):
    if frac <= 0.0 + 1e-9: return 0, 0
    elif frac <= 0.075 + 1e-9: return 1, 0
    else: return 1, 1


def bin_to_gray(x): 
    return x ^ (x >> 1)

def repair_individual(g):
    def repair_enable_connectivity(g):
        # 把 MST 边强制打开，作为连通骨架
        en = g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1]
        for idx in mst_indices:
            en[idx] = 1
        g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1] = en
    
    def repair_monotone_diam_codes(g):
        # 1) 读成索引（0..9），先 Gray→bin，再截断
        idxs = np.zeros(NUM_PIPES, dtype=int)
        for i in range(NUM_PIPES):
            off = IDX_DIAM_START + 4*i
            code = (g[off] | (g[off+1]<<1) | (g[off+2]<<2) | (g[off+3]<<3))
            idxs[i] = min(gray_to_bin(int(code)), 9)

        # 2) 抬下游到不小于上游最大
        changed = True
        while changed:
            changed = False
            for (u,v), di in edge_to_pipe_idx.items():
                if not g[IDX_PIPE_EN_START + di]: 
                    continue
                # 所有 w->u 的上游
                for (w,uu), ui in edge_to_pipe_idx.items():
                    if uu == u and g[IDX_PIPE_EN_START + ui]:
                        if idxs[ui] > idxs[di]:
                            idxs[di] = idxs[ui]; 
                            changed = True

        # 3) 回写 4 位 Gray
        for i in range(NUM_PIPES):
            code = bin_to_gray(int(idxs[i]))
            off = IDX_DIAM_START + 4*i
            g[off+0] = (code >> 0) & 1 # 小端
            g[off+1] = (code >> 1) & 1
            g[off+2] = (code >> 2) & 1
            g[off+3] = (code >> 3) & 1
    def repair_outlet_bits_M(g):
        '''
        # 读启用边 → 子图
        enabled = g[IDX_PIPE_EN_START:IDX_PIPE_EN_END+1].astype(np.uint8)
        H = nx.DiGraph(); H.add_nodes_from(G_full.nodes())
        H.add_edges_from([(u, v) for (u, v), idx in edge_to_pipe_idx.items() if enabled[idx]])
        # 候选末端（排除子流域接口）
        sinks = [n for n in H.nodes() if H.out_degree(n) == 0 and (subcatch_outlets is None or n not in subcatch_outlets)]
        chosen = []
        if len(sinks) >= 2:
            # 也可用 random.sample(sinks, 2)
            chosen = sinks[:2]
        else:
            # 末端不够则从 node_list 顺序补齐
            chosen = sinks.copy()
            for n in node_list:
                if n not in chosen:
                    chosen.append(n)
                if len(chosen) == 2: break

        # 清零出口段 → 恢复恰好两个 1
        outM = np.zeros(M, dtype=np.uint8)
        for n in chosen[:2]:
            outM[node_to_idx[n]] = 1
        '''
        outM = OUTLET_FIXED_MASK.copy()  # —— 关键：写入 M
        g[IDX_OUTLET_START:IDX_OUTLET_END+1] = outM
        return g      
    
    def repair_lid_bits(g):
        # 逐子汇：读取 2bit→frac，若 sum>0.10 则截断并回写
        for i in range(NUM_SUBCATCH):
            en_bc = int(g[IDX_BC_START + 2*i + 0])
            ab_bc = int(g[IDX_BC_START + 2*i + 1])
            en_pp = int(g[IDX_PP_START + 2*i + 0])
            ab_pp = int(g[IDX_PP_START + 2*i + 1])

            bc = _bits_to_frac(en_bc, ab_bc)
            pp = _bits_to_frac(en_pp, ab_pp)
            s  = bc + pp
            if s > 0.10 + 1e-9:
                # 简单策略：优先保留较大的那一项，把另一项减到恰好 0.10
                if bc >= pp:
                    bc = min(bc, 0.10); pp = max(0.0, 0.10 - bc)
                else:
                    pp = min(pp, 0.10); bc = max(0.0, 0.10 - pp)
                # 回写为 2bit
                en_bc, ab_bc = _frac_to_bits(bc)
                en_pp, ab_pp = _frac_to_bits(pp)
                g[IDX_BC_START + 2*i + 0] = en_bc
                g[IDX_BC_START + 2*i + 1] = ab_bc
                g[IDX_PP_START + 2*i + 0] = en_pp
                g[IDX_PP_START + 2*i + 1] = ab_pp  
    
    # 先确保基本连通
    repair_enable_connectivity(g)
    # 修复管径单调
    repair_monotone_diam_codes(g)
    # 修复出口（恰好 2 个末端）
    repair_outlet_bits_M(g)
    # 修复 LID（BC+PP ≤ 0.10）
    repair_lid_bits(g)
    return g


def generate_inp_file(template_inp_path, output_inp_path,
                      pipe_enabled, pipe_diameters, conduits_info,
                      outlet_nodes,
                      lid_vars,
                      loops,          # numpy array from Loops.txt: each row lists pipe indices in that subcatchment
                      subcatch_info, # list or array of areas for S1…S18
                      coords_info, # list of coordinates for each node
                      rain_gage="RG1"):
    """
    Rewrite only these blocks:
      [OUTFALLS]
      [CONDUITS]
      [XSECTIONS]
      [SUBCATCHMENTS]
      [SUBAREAS]
      [LID_CONTROLS]
      [LID_USAGE]
    Leave everything else intact.
    """
    lines = open(template_inp_path).read().splitlines(True)

    def replace_block(lines, section, new_block):
        tgt = f'[{section.upper()}]'
        # find section header
        for i,l in enumerate(lines):
            if l.strip().upper() == tgt:
                start = i
                break
        else:
            raise RuntimeError(f"Section [{section}] not found")
        # find end of block
        end = start+1
        while end < len(lines) and not lines[end].strip().startswith('['):
            end += 1
        return lines[:start+1] + new_block + lines[end:]

    
    # --- Build the directed graph to find outlets ---
    #outlets = find_outlets(pipe_enabled)

    # 1) OUTFALLS
    outfalls = [";;Name           Elevation  Type       Stage Data       Gated    Route To\n",        
                ";;-------------- ---------- ---------- ---------------- -------- ----------------\n"]
    for sel in system_outfall_node:
        elev = elev_map.get(sel, 0.0)
        outfalls.append(
            f"{sel:<15d}"
            f"{elev:>8.2f}   "
            f"{'FREE':<9s}"
            f"{'':<11s}"
            f"{'NO':<8s}"
            f"\n"
        )
    lines = replace_block(lines, "OUTFALLS", outfalls)

    # 2) CONDUITS
    conduits = [";;Name           From Node        To Node          Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow\n",
                ";;-------------- ---------------- ---------------- ---------- ---------- ---------- ---------- ---------- ----------\n"] 
    # 我们要确保：每个 outfall_selection 中的 node 只被一条管道指向；而且绝不出现在 from_node
    used_outfalls = set() 
    #emitted = []              # <— keep track of which pipe_ids were actually written
    for row in conduits_info:
        pipe_id = int(row[0])
        frm = int(row[1]); to = int(row[2])
        # 1) 只保留那些被允许的管道
        idx = edge_to_pipe_idx.get((frm, to))
        if idx is None or not pipe_enabled[idx]:
            continue
        # 2) 如果这个管道的 to_node 恰好是系统 outfall
        #    且我们已经用过一次了，就跳过它
        if to in outlet_nodes:
            if to in used_outfalls:
                # 禁用掉第二条以上接到同一节点的管线
                continue
            # 第一次遇到，记录
            used_outfalls.add(to)

        # 3) 同时，确保这个 outfall 也不会被当作 from_node
        if frm in outlet_nodes:
            # 如果它是 outfall，就永远不能做管线的起点
            continue

        L = float(row[3]); R = float(row[4])
        conduits.append(
            f"{pipe_id}\t{frm}\t{to}\t{L:>8.2f}\t{R:>8.4f}\t0\t0\t0\t0\n"
        )
        #emitted.append((j, pipe_id))
    lines = replace_block(lines, "CONDUITS", conduits)

    # 3) XSECTIONS (holds the diameters)
    xsecs = [";;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels    Culvert\n",   
             ";;-------------- ------------ ---------------- ---------- ---------- ---------- ---------- ----------\n"]
    used_outfalls = set()
    for row in conduits_info:
        pipe_id = int(row[0])
        frm     = int(row[1])
        to      = int(row[2])
        # 找全局 idx，看它是不是 enabled
        idx = edge_to_pipe_idx.get((frm, to))
        if idx is None or not pipe_enabled[idx]:
            continue

        # 同样的：一个 outfall 只能被用一次
        if to in system_outfall_node:
            if to in used_outfalls:
                continue
            used_outfalls.add(to)

        # 出水口不能出度
        if frm in system_outfall_node:
            continue
        d = pipe_diameters[idx]
        
        if d is None or d <= 0:
            # print when d is None or <= 0
            print(f"Warning: Pipe {pipe_id} has invalid diameter {d}. Using default 0.25 m.")
            d = 0.25
        xsecs.append(f"{pipe_id:<15d}\tCIRCULAR\t{d:.2f}\t0\t0\t0\t1\t\t\n")
    lines = replace_block(lines, "XSECTIONS", xsecs)
    # print xsecs
    #print("DEBUG: XSECTIONS block:", xsecs)
    


    # 4) SUBCATCHMENTS
    subcat = [
        ";;Name           Rain Gage        Outlet           Area     %Imperv  Width    %Slope   CurbLen  SnowPack        \n",
        ";;-------------- ---------------- ---------------- -------- -------- -------- -------- -------- ----------------\n"
    ]
    # subcatch_info.shape[0] should be 18 (one line per catchment)
    
    for i in range(n_sc):
        sid     = i+67         # 子流域编号，如 1→S1
        outlet  = int(subcatchout_info[i, 1])         # 模板期望的出口节点号
        row      = subcatch_info[i]
        rg      = int(row[1])         # 雨量计编号
        area    = row[3]              # 面积
        imperv  = row[4]              # 不透水率（%）
        width   = row[5]              # 宽度
        slope   = row[6]              # 坡度（%）
        curb    = row[7]              # 路缘长度
        #snow    = row[8]              # 雪堆积（mm）
        # SnowPack we fill with 0
        subcat.append(
            f"{sid:<14d}"   # 左对齐，占14列
            f"{rg:<16d}"    # 左对齐，占16列
            f"{outlet:<16d}"# 左对齐，占16列
            f"{area:>8.2f}" # 右对齐，占8列，保留2位小数
            f"{imperv:>8.2f}"
            f"{width:>8.2f}"
            f"{slope:>8.2f}"
            f"{curb:>8.2f}"
            f"\n"            # (空) SnowPack 列
        )
    lines = replace_block(lines, "SUBCATCHMENTS", subcat)

    # 5) SUBAREAS
    subareas = [";;Subcatchment   N-Imperv   N-Perv     S-Imperv   S-Perv     PctZero    RouteTo    PctRouted \n",
                ";;-------------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n"
    ]
    
    # subcatchment: 67 to 84
    for sid in range(67, 85):
        subareas.append(
            f"{sid}\t0.024\t0.15\t2.1\t6.51\t30\tOUTLET\n"
        )
    lines = replace_block(lines, "SUBAREAS", subareas)

    # 6) INFILTRATION
    # fixed parameters for every subcatchment
    max_rate  = 103.81
    min_rate  =  11.44
    decay     =   2.75
    dry_time  =   7
    max_infil =   0
    infiltration = [";;Subcatchment   MaxRate    MinRate    Decay      DryTime    MaxInfil \n", 
                    ";;-------------- ---------- ---------- ---------- ---------- ----------\n"]
    for sid in range(67, 85):
        # sid is the new ID (67…84) you want in column 1
        infiltration.append(
            f"{int(sid):<14d} "           # Subcatchment ID
            f"{max_rate:>10.2f} "         # MaxRate
            f"{min_rate:>10.2f} "         # MinRate
            f"{decay:>10.2f} "            # Decay
            f"{dry_time:>10d} "           # DryTime
            f"{max_infil:>10d}\n"         # MaxInfil
        )
    lines = replace_block(lines, "INFILTRATION", infiltration)
    
    # 7) LID_CONTROLS
    lidc = [
    ";;Name           Type/Layer Parameters\n",
    ";;-------------- ---------- ----------\n",
    # First unit (1 = BC)
    "1                BC\n",
    "1                SURFACE    450        0.05       0.1        0.5        5         \n",
    "1                SOIL       900        0.5        0.15       0.08       50         10         80        \n",
    "1                STORAGE    300        0.67       500         0          \n",
    "1                DRAIN      2.5        2.5        150        6          0          0         \n",
    "\n",
    # Second unit (2 = PP)
    "2                PP\n",
    "2                SURFACE    0.0        0.0        0.012      0.5        5         \n",
    "2                PAVEMENT   100        0.15       0          500        0          0          0         \n",
    "2                STORAGE    300        0.4        500         0          \n",
    "2                DRAIN      2.5        0.5        100        6          0          0         \n",
]
    
    lines = replace_block(lines, "LID_CONTROLS", lidc)

    # 8) LID_USAGE
    lidu = [";;Subcatchment   LID Process      Number  Area       Width      InitSat    FromImp    ToPerv     RptFile                  DrainTo\n",         
            ";;-------------- ---------------- ------- ---------- ---------- ---------- ---------- ---------- ------------------------ ----------------\n"]
    start_id = 67
    width = 10
    unit_area_bc  = 200     # 每个 BC 单元理论最大面积
    unit_area_pp  = 150 
    fromImp_BC = 57.14286
    fromImp_PP = 42.85714
    bc_fracs = lid_vars[0::2]   # length 18, fraction of each catchment’s area used by BC
    pp_fracs = lid_vars[1::2]   # length 18, fraction for PP
    for i, (bc_frac, pp_frac, A) in enumerate(zip(bc_fracs, pp_fracs, catchment_areas)):
        
        sid = start_id + i
        # --- BC 单元 ---
        # 这里我们假设每个子流域最多只有一个 BC 单元 (Number=1),
        # 如果你想动态设置 num_bc，可以从 lv 里读 lv["num_bc"]
        bc_area = bc_frac * A
        if bc_area > 0:
            num_bc       = int(np.ceil(bc_area / unit_area_bc)) if bc_area > 0 else 0
            #area_per_bc  = bc_area / num_bc
            from_imp_bc  = fromImp_BC if num_bc > 0 else 0.0
            # 第七列 FromImp = 单元面积 / 子流域面积 * 100 (%)
            lidu.append(
                f"{sid:<14d}"      # Subcatchment ID, 左对齐 14 列
                f"{1:<16d}"        # LID Process: 对 BC 赋码 1
                f"{num_bc:<9d}"    # Number
                f"{unit_area_bc:>10.2f}"# Area
                f"{width:>11d}"    # Width
                f"{0:>11d}"        # InitSat
                f"{from_imp_bc:>11.5f}"# FromImp
                f"{1:>11d}"        # ToPerv
                "\n"
            )
        # --- PP 单元 ---
        pp_area = pp_frac * A
        if pp_area > 0:
            num_pp       = int(np.ceil(pp_area / unit_area_pp)) if pp_area > 0 else 0
            #area_per_pp  = pp_area / num_pp
            from_imp_pp  = fromImp_PP if num_pp > 0 else 0
            lidu.append(
                f"{sid:<14d}"
                f"{2:<16d}"        # LID Process: 对 PP 赋码 2
                f"{num_pp:<9d}"
                f"{unit_area_pp:>10.2f}"
                f"{width:>11d}"
                f"{0:>11d}"
                f"{from_imp_pp:>11.5f}"
                f"{1:>11d}"
                "\n"
            )
    #print("DEBUG lid_vars:", lid_vars)      
    # 把 LID_USAGE 区块写回去
    lines = replace_block(lines, "LID_USAGE", lidu)
  
    # 9) Junctions
    
    junc_block = [
        ";;Name           Elevation  MaxDepth   InitDepth  SurDepth   Aponded   \n",
        ";;-------------- ---------- ---------- ---------- ---------- ----------\n"
    ]
    for row in junctions_info:
        node = int(row[0])
        if node in system_outfall_node:
            # 如果是出水口节点，跳过
            continue   # skip outlet nodes here
        elev     = row[1]
        maxD     = row[2]
        initD    = row[3]
        surD     = row[4]
        aponded  = row[5]
        junc_block.append(
            f"{node:<7d}{elev:>11.4f}{maxD:>11.4f}"
            f"{initD:>11.4f}{surD:>11.4f}{aponded:>11.4f}\n"
        )
    lines = replace_block(lines, "JUNCTIONS", junc_block)
    
    # 10) Coordinates
    #co_sc = coords_info.shape[0] # subcatchment ID
    coord = [";;Node           X-Coord            Y-Coord\n",
             ";;-------------- ------------------ ------------------\n"]           
    for node, x, y in coords_info:
        coord.append(
            f"{int(node):<16d}"   # Node ID, left‑aligned in 16 chars
            f"{x:>12.1f}"         # X‑Coord, right‑aligned, 1 decimal
            f"{y:>13.1f}"         # Y‑Coord, right‑aligned, 1 decimal
            "\n"
        )
    lines = replace_block(lines, "COORDINATES", coord)
    
    # 11) Polygons
    poly = [";;Subcatchment   X-Coord            Y-Coord\n",           
            ";;-------------- ------------------ ------------------\n",
            "67       1962.838       7186.937\n",
            "68       4902.872       6897.523\n",
            "69       3799.268        6604.73\n",
            "70       3101.914       6355.856\n",
            "71       2608.671       5902.027\n",
            "72       4062.853       5361.486\n",
            "73       5137.318       4888.795\n",
            "74       5865.218       4698.198\n",
            "75       6804.407       4385.135\n",
            "76       7749.227       4034.347\n",
            "77       8664.205       3744.369\n",
            "78       9278.716       3402.027\n",
            "79        10071.79        2957.489\n",
            "80       11300.676        2972.973\n",
            "81        9374.87       4591.377\n",
            "82       10647.393        4552.767\n",
            "83       9373.825       5799.549\n",
            "84       9734.889       5242.117\n"]

    lines = replace_block(lines, "POLYGONS", poly)
    
    # write out
    with open(output_inp_path, 'w') as f:
        f.writelines(lines)

    return output_inp_path

def get_system_qmax_from_rpt(rpt_path):
    """从已有的 rpt 文件里提取 System 最大流量，不跑 SWMM。"""
    with open(rpt_path) as f:
        for raw in f:
            line = raw.strip()
            # 忽略大小写，只要行以 System 开头就解析第 4 列
            if re.match(r'(?i)^system\b', line):
                parts = re.split(r'\s+', line)
                if len(parts) >= 4:
                    try:
                        return float(parts[3])
                    except ValueError:
                        pass
    raise RuntimeError(f"未在 {rpt_path} 找到 System 最大流量。")

def compute_LCC(inp_path,cost_per_diameter,lid_vars,cost_lid_bc=200,cost_lid_pp=150,om_lid_bc=0.08,om_lid_pp=0.04,discount_rate=0.02,project_lifetime=30):
    """
    1) 读 [XSECTIONS]：第三列为直径 → dict name->diam
    2) 读 [CONDUITS]：第四列为 Length，第一列为 Link Name
       cost_gray = sum(length_i * cost_per_diameter[diam_i])
    3) cost_bc = cost_lid_bc * ΣBC_area
       cost_pp = cost_lid_pp * ΣPP_area
    4) om_cost_bc = cost_bc * om_lid_bc * PV_factor
       om_cost_pp = cost_pp * om_lid_pp * PV_factor
    PV_factor = (1 - (1+discount_rate)**(-project_lifetime)) / discount_rate
    """
    lcc = 0.0

    # --- 1) 解析 XSECTIONS ---
    xs_diam = {}
    with open(inp_path, 'r') as f:
        in_xs = False
        for line in f:
            line = line.strip()
            if not in_xs:
                if line.upper().startswith("[XSECTIONS]"):
                    in_xs = True
                    continue
            else:
                if line.startswith("["):
                    break
                if line and not line.startswith(";;"):
                    parts = line.split()
                    name = parts[0]
                    # parts[2] 是第三列直径
                    diam = float(parts[2])
                    xs_diam[name] = diam

    # --- 2) 解析 CONDUITS 计算灰设施成本 ---
    cost_gray = 0.0
    with open(inp_path, 'r') as f:
        in_cd = False
        for line in f:
            line = line.strip()
            if not in_cd:
                if line.upper().startswith("[CONDUITS]"):
                    in_cd = True
                    continue
            else:
                if line.startswith("["):
                    break
                if line and not line.startswith(";;"):
                    parts = line.split()
                    link = parts[0]
                    length = float(parts[3])   # 4th colunm, Length
                    diam = xs_diam.get(link)
                    if diam is None:
                        # if nothing in XSECTIONS, skip 
                        continue
                    unit_cost = cost_per_diameter.get(diam, 0.0)
                    cost_gray += length * unit_cost

    # --- 3) LID 造价（一次性投资） ---
    bc_fracs = lid_vars[0::2]
    pp_fracs = lid_vars[1::2]
    total_bc = (bc_fracs * catchment_areas).sum()
    total_pp = (pp_fracs * catchment_areas).sum()
    cost_bc = cost_lid_bc * total_bc
    cost_pp = cost_lid_pp * total_pp

    # --- 4) LID 维护运维成本折现 ---
    pv_factor = (1 - (1 + discount_rate) ** (-project_lifetime)) / discount_rate
    om_cost_bc = cost_bc * om_lid_bc * pv_factor
    om_cost_pp = cost_pp * om_lid_pp * pv_factor

    lcc = cost_gray + cost_bc + cost_pp + om_cost_bc + om_cost_pp
    return lcc
    

def calculate_fitness(individual, LCC_MIN, LCC_MAX, Q_MAX_MIN, Q_MAX_MAX, inp_path, rpt_path):
    """
    计算适应度函数，结合 LCC 和最大径流 Q_max。
    """
    # 1) 解码
    pipe_enabled, pipe_diameters, outlet_nodes, lid_vars = decode_genome_M(individual, node_list)

    
    # 可行性预检失败 → 尝试修复一次，再 decode & 复检
    if (not check_lid_constraints(lid_vars, catchment_areas, 0.10)) or (not check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes, edge_to_pipe_idx, G_full, subcatch_outlets)):
        individual = repair_individual(individual.copy())
        pipe_enabled, pipe_diameters, outlet_nodes, lid_vars = decode_genome_M(individual, node_list)
        # 再检一次，不通过再罚
        if (not check_lid_constraints(lid_vars, catchment_areas, 0.10)) or \
        (not check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes, edge_to_pipe_idx, G_full, subcatch_outlets)):
            return -1.0
    
    try:
        ok_struct = check_structure_validity(
            pipe_enabled, pipe_diameters, outlet_nodes,
            edge_to_pipe_idx=edge_to_pipe_idx, G_full=G_full,
            subcatch_outlets=subcatch_outlets
        )
    except Exception:
        ok_struct = False
    if not ok_struct:
        return -2.0

    # 3) 生成 INP（注意传入 outlet_nodes）
    tmodify = time.time()
    generate_inp_file(
        template_inp_path, inp_path,
        pipe_enabled, pipe_diameters, conduits_info,
        outlet_nodes,                      # ← 新增
        lid_vars,
        loops, subcatch_info, coords_info
    )
    # print(f"[FIT] modify_inp_file: {time.time()-tmodify:.2f}s")

    # 4) 跑 SWMM + 取 LCC/Qmax
    tlcc = time.time()
    lcc = compute_LCC(inp_path, cost_per_diameter, lid_vars,
                  cost_lid_bc=200, cost_lid_pp=150,
                  om_lid_bc=0.08, om_lid_pp=0.04,
                  discount_rate=0.02, project_lifetime=30)
    # qmax from swmm and rpt
    subprocess.run(
        [swmm_exe_path, inp_path, rpt_path],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    # print(f"[FIT] calculate_lcc: {time.time()-tlcc:.2f}s")

    q_max = get_system_qmax_from_rpt(rpt_path)

    # 5) 归一化 + 适应度（越大越好）
    lcc_norm = (lcc  - LCC_MIN)   / (LCC_MAX   - LCC_MIN   + EPS)
    q_norm   = (q_max - Q_MAX_MIN) / (Q_MAX_MAX - Q_MAX_MIN + EPS)
    fitness  = w1*(1 - lcc_norm) + w2*(1 - q_norm)
    return float(fitness)

# Selection: Tournament Selection
def tournament_selection(population, fitness_values, k=3):
    """
    Selects a parent using tournament selection.
    - k: Tournament size (default: 3)
    """
    selected = np.random.choice(len(population), k, replace=False)
    best_index = selected[np.argmax(fitness_values[selected])]
    return population[best_index]

# Crossover: One-Point Crossover
def crossover(parent1, parent2):
    """
    Performs single-point crossover.
    """
    if random.random() < crossover_rate:
        point = random.randint(1, TOTAL_BITS - 1)  # Choose crossover point
        child1 = np.concatenate((parent1[:point], parent2[point:]))
        child2 = np.concatenate((parent2[:point], parent1[point:]))
    else:
        child1, child2 = parent1.copy(), parent2.copy()  # No crossover
    return child1, child2

# Mutation: Bit Flip Mutation
def mutate(individual):
    """
    Performs mutation by flipping random bits.
    """
    for i in range(TOTAL_BITS):
        if random.random() < mutation_rate:
            individual[i] = 1 - individual[i]  # Flip bit
    return individual


def genetic_algorithm(inp_path, rpt_path):
    """
    Runs Genetic Algorithm to optimize infrastructure design.
    """
    global fitness_trend, lcc_trend, qmax_trend  # 记录趋势
    time_trend = []
    fitness_trend = []  # 清空趋势数据
    lcc_trend = []       # 每代最佳 LCC
    qmax_trend = []      # 每代最佳 Q_max
    norm_obj_trend = []
    lcc_norm_trend = []
    q_norm_trend = []

    
    start_time = time.time() # start timing
    
    # Initialize population
    print("begin to initialize population...")
    t0 = time.time()
    population = init_population_feasible(population_size)
    assert population.ndim == 2 and population.shape[1] == TOTAL_BITS, \
        f"population.shape[1]={population.shape[1]} != TOTAL_BITS={TOTAL_BITS}"
    print(f"[POP OK] shape={population.shape}, TOTAL_BITS={TOTAL_BITS}")
    print(f"[GA] 初始化种群耗时: {time.time()-t0:.2f}s")
    
    # 设置新 `.inp` 文件路径（为 GA 生成的新文件）
    #new_inp = "./modifiedswmm.inp"
    #new_rpt = "./modifiedswmm.rpt"

    # Compute initial LCC and Max Outflow range for normalization
    print("begin to calculate LCC and Max Outflow...")
    
    '''
    # used the designated range instead of calculated range, that's why this part is commented
    lcc_values = np.array([calculate_lcc(individual) for individual in population])
    outflow_values = np.array([
        run_swmm_simulation(new_inp) if (individual.sum() > 0) else 0.0  # avoid empty individual
        for individual in population
    ])
    LCC_MIN, LCC_MAX = np.min(lcc_values), np.max(lcc_values)
    Q_MAX_MIN, Q_MAX_MAX = np.min(outflow_values), np.max(outflow_values)
    '''
    
    #LCC_MIN, LCC_MAX = 1e5, 3e6
    #Q_MAX_MIN, Q_MAX_MAX = 0, 3e4
    

    print("begin to evaluate population...")   
    t1 = time.time()
    # Evaluate initial population
    fitness_values = np.array([
        calculate_fitness(individual, LCC_MIN, LCC_MAX, Q_MAX_MIN, Q_MAX_MAX, inp_path, rpt_path)
        for individual in population
    ])
    print(f"[GA] 初始人群评估耗时: {time.time()-t1:.2f}s")

    for generation in range(num_generations):
        print(f"Generation {generation+1}/{num_generations}...")
        # set up a limit for generation, if needed
        #if generation >= 10:  # enforced termination after 10 generations
        #    break
        
        new_population = []
        t2 = time.time()

        '''
        # **Elitism: Retain best individuals** (temporarily disabled)
        elite_indices = np.argsort(fitness_values)[-elite_size:]
        elite_individuals = population[elite_indices]

        # Generate new population
        while len(new_population) < population_size - elite_size:
            print("tournament selection...")
            parent1 = tournament_selection(population, fitness_values)
            parent2 = tournament_selection(population, fitness_values)
            child1, child2 = crossover(parent1, parent2)
            new_population.append(mutate(child1))
            new_population.append(mutate(child2))
        
            
        # Ensure new population size matches original
        new_population = np.array(new_population[:population_size - elite_size])
        new_population = np.vstack((elite_individuals, new_population))  # Keep elite
        '''
        while len(new_population) < population_size:
            parent1 = tournament_selection(population, fitness_values)
            parent2 = tournament_selection(population, fitness_values)
            child1, child2 = crossover(parent1, parent2)
            child1 = mutate(child1); child1 = repair_individual(child1)
            child2 = mutate(child2); child2 = repair_individual(child2)
            new_population.append(child1)
            new_population.append(child2)   

        # Ensure new population size matches original
        population = np.array(new_population[:population_size])


        # Update population
        print("update population...")
        print(f"[GA] 选择/交叉/变异耗时(第{generation+1}代): {time.time()-t2:.2f}s")

        # Recalculate fitness
        print("recalculate fitness...")
        t3 = time.time()
        fitness_values = np.array([
            calculate_fitness(individual, LCC_MIN, LCC_MAX, Q_MAX_MIN, Q_MAX_MAX, inp_path, rpt_path)
            for individual in population
        ])
        print(f"[GA] 评估整群耗时(第{generation+1}代): {time.time()-t3:.2f}s")
        
        # Check for convergence acording to fitness values, if needed
        #if np.std(fitness_values) < 1e-6:
        #    print("⚠️ Fitness 收敛过快，跳出 GA 提前终止。")
        #    break
        
        # print max outflow each generation
        t4 = time.time()
        #max_outflow = max([run_swmm_simulation(new_inp) for individual in population])
        #print(f"Generation {generation+1}: Max Outflow = {max_outflow:.4f}")
        print(f"[GA] 扫描最大出流耗时(第{generation+1}代): {time.time()-t4:.2f}s")
        
        # 记录最优（针对最优个体重写INP并跑 SWMM）
        best_index = np.argmax(fitness_values)
        best_individual = population[best_index]
        # Print best fitness, LCC and outflow per generation
        best_fitness = np.max(fitness_values)
        if best_fitness <= -1 + 1e-9:
            # 打印约束不满足的具体原因
            debug_constraints_for_individual(
                best_individual,
                node_list=node_list,
                catchment_areas=catchment_areas,
                edge_to_pipe_idx=edge_to_pipe_idx,
                G_full=G_full,
                subcatch_outlets=subcatch_outlets,
                require_single_component=True,      # 或 False 试试放宽
                require_reachability=False
            )
        fitness_trend.append(best_fitness)
        print(f"Generation {generation+1}: Best Fitness = {best_fitness:.4f}")
        elapsed = time.time() - start_time
        time_trend.append(elapsed)
        
        t5 = time.time()
        pipe_enabled, pipe_diameters, outlet_nodes, lid_vars = decode_genome_M(best_individual, node_list)
        generate_inp_file(template_inp_path, inp_path, pipe_enabled, pipe_diameters, conduits_info, outlet_nodes, lid_vars, loops, subcatch_info, coords_info)
        best_lcc = compute_LCC(inp_path, cost_per_diameter, lid_vars,
                      cost_lid_bc=200, cost_lid_pp=150,
                      om_lid_bc=0.08, om_lid_pp=0.04,
                      discount_rate=0.02, project_lifetime=30)
        subprocess.run([swmm_exe_path, inp_path, rpt_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        best_qmax = get_system_qmax_from_rpt(rpt_path)
        lcc_trend.append(best_lcc)
        qmax_trend.append(best_qmax)      
        
        # 与 calculate_fitness 里同样的方法做归一化
        lcc_norm   = (best_lcc  - LCC_MIN)   / (LCC_MAX   - LCC_MIN   + EPS)
        q_norm     = (best_qmax - Q_MAX_MIN) / (Q_MAX_MAX - Q_MAX_MIN + EPS)
        #lcc_norm   = min(max(lcc_norm, 0), 1)
        #q_norm     = min(max(q_norm,   0), 1)
        
        # 计算“归一化目标值”（这里跟 fitness 定义严格一致）
        J_norm     = w1 * lcc_norm + w2 * q_norm
        norm_obj_trend.append(J_norm)
        lcc_norm_trend.append(lcc_norm)
        q_norm_trend.append(q_norm)
        print(f"[GA] 记录最佳个体指标耗时(第{generation+1}代): {time.time()-t5:.2f}s")
        # print norm_trend of each obj
        print(f"Generation {generation+1}: LCC Norm = {lcc_norm:.4f}, Qmax Norm = {q_norm:.4f}, Objective = {J_norm:.4f}, best LCC = {best_lcc:.2f}, best Qmax = {best_qmax:.2f}")
        
    # Best solution found
    print("\nBest Individual Found:", best_individual)
    print("Best LCC:", best_lcc)
    #print("Best Max Outflow:", run_swmm_simulation(new_inp))
    print("Best Max Outflow:", best_qmax)
    print("Best Fitness:", fitness_values[best_index])
    
    end_time = time.time()
    duration = end_time - start_time
    print(f"🕒Time duration: {duration:.2f} seconds" )    
        
    # plot convergence trends
    plt.subplot(1, 3, 1)
    plt.plot(range(1, len(norm_obj_trend)+1), norm_obj_trend, marker="d", linestyle="-")
    plt.ylim(0, 1)
    plt.xlabel("Generation")
    plt.ylabel("Normalized Objective")
    plt.title("Objective Convergence")
    plt.grid(True)
    
    # plot LCC
    plt.subplot(1, 3, 2)
    plt.plot(range(1, len(lcc_norm_trend) + 1), lcc_norm_trend, marker="s", linestyle="--", color="green", label="LCC")
    plt.ylim(bottom=0.5 * min(lcc_norm_trend))  # y轴下限从合理最小值开始
    plt.xlabel("Generation")
    plt.ylabel("LCC")
    plt.title("LCC Convergence")
    plt.legend()
    plt.grid(True)
    
    # plot Qmax
    plt.subplot(1, 3, 3)
    plt.plot(range(1, len(q_norm_trend) + 1), q_norm_trend, marker="^", linestyle="--", color="orange", label="Qmax")
    plt.ylim(bottom=0.5 * min(q_norm_trend))  # y轴下限从合理最小值开始
    plt.xlabel("Generation")
    plt.ylabel("Max Outflow (CMS)")
    plt.title("Qmax Convergence")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(6,4))
    plt.plot(time_trend, norm_obj_trend, marker='o', linestyle='-')
    plt.xlabel("Elapsed Time (s)")
    plt.ylabel("Best Objective Value")
    plt.title("Objective Convergence vs Time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    template_inp_path = "./AhvazNull.inp"
    swmm_exe_path = "./swmm5.exe"
    inp_path = "./modifiedswmm.inp"
    rpt_path = "./modifiedswmm.rpt"
    #base_inp = "./LIMatF2.inp"
    #out_file = "./LIMatF2.out"
    # 生成 1:9 → 9:1 的 9 组权重
    #weights = [(i/10.0, 1.0 - i/10.0) for i in range(1, 10)]
    total_t0 = time.time()  
    best_solution = genetic_algorithm(inp_path, rpt_path)
    #run_for_weights_and_plot(base_inp, weights)
    #print(f"[TOTAL] 全流程耗时: {time.time()-total_t0:.2f}s")



# ga_swmm_binary_329.py
# -------------------------------------------------------------
# Standalone Genetic Algorithm (binary genome, total 329 bits)
# Genome layout (0-indexed, inclusive ranges):
#   Pipe enable bits:      [  0 ..  82]  -> 83 bits (each 1/0: keep/remove a pipe)
#   Outlet select bits:    [ 83 ..  90]  -> 8 bits (each 1/0: enable a potential outfall node)
#   Pipe diameter bits:    [ 91 .. 256]  -> 166 bits = 83 pipes × 2 bits per pipe
#   LID BC bits:           [257 .. 292]  -> 36 bits = 18 subcatch × (enable 1 bit + area 1 bit)
#   LID PP bits:           [293 .. 328]  -> 36 bits = 18 subcatch × (enable 1 bit + area 1 bit)
# Total = 329 bits
#
# Notes:
# - Independent script; no external project files required.
# - Uses a surrogate hydraulic/cost model by default (fast, portable).
# - You can later hook it to SWMM by replacing the surrogate in `evaluate()`.
# - All parameters (weights, economics, GA hyperparams) are editable below.
# -------------------------------------------------------------












































