import numpy as np
import pandas as pd
import random
import torch
import gpytorch
import subprocess
import re
import os
import networkx as nx
from pyswmm import Output, Simulation, Nodes
import time
import matplotlib.pyplot as plt
import random
from typing import List, Dict, Tuple, Optional, Sequence   
from collections import defaultdict
from pyDOE2 import lhs

# Number of initial samples and BO iterations
N_initial = 2000
max_iter = 100         

# Normalization bounds
LCC_MIN, LCC_MAX = 3e4, 3e6
Q_MAX_MIN, Q_MAX_MAX = 0, 3e4
EPS = 1e-6

# Objective weights (minimize w1*LCC_norm + w2*Qmax_norm)
w1 = 0.5  # weight for LCC
w2 = 0.5  # weight for max outflow

# Load conduits
conduits_info = np.loadtxt("[CONDUITS].txt")

# Load commercial diameters

commercial_df = pd.read_csv("Commercial.txt", sep=r"\s+", header=None)
commercial_df.columns = ["Diameter_m", "UnitCost"]
commercial_diameters = commercial_df["Diameter_m"].tolist()
# cost dictionary from file (if needed)
cost_per_diameter = dict(zip(commercial_df["Diameter_m"], commercial_df["UnitCost"]))

        

# 离散可选权重（管径）列表，严格递增
#weight_list = commercial_diameters.copy()
#weight_to_state = {w: commercial_diameters.index(w) + 1 for w in commercial_diameters}
weight_list = commercial_diameters.copy()
weight_to_state = { w: weight_list.index(w) + 1 for w in weight_list }
cost_per_diameter = {float(w): float(cost_per_diameter[float(w)]) for w in weight_list}

loops = np.loadtxt("Loops.txt", dtype=int)
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
M = len(node_list)      # 出口槽数: 48 
LID_LEVELS = (0.0, 0.025, 0.05, 0.10)


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
    
D = N_pipes + M + 36 # 设计向量总维度 (dynamic)
#d_cont = 36
#d_disc = 65 + 18 + M
d_cont = N_pipes 
d_disc = M + 36

def repair_monotone_global(pipe_enabled: np.ndarray,
                           pipe_diameters: List[Optional[float]],
                           weight_list: List[float],
                           edge_to_pipe_idx: Dict[Tuple[int, int], int],
                           n_iter: int = 10) -> None:
    """全网‘上游≤下游’修复：把直径先映射为离散 code，再多轮抬高下游。原地修改 pipe_diameters。"""
    L = len(weight_list)
    code_of = {w: i for i, w in enumerate(weight_list)}
    codes = np.full(len(pipe_enabled), -1, dtype=int)
    for i, (en, d) in enumerate(zip(pipe_enabled, pipe_diameters)):
        if en and (d in code_of):
            codes[i] = code_of[d]
    for _ in range(n_iter):
        changed = False
        # 注意：这里假设 edge_to_pipe_idx 的 key 里包含所有有向边 (u,v)
        for (u, v), idx_d in edge_to_pipe_idx.items():
            if not pipe_enabled[idx_d] or codes[idx_d] < 0:
                continue
            # 所有“直接上游”的边 (w,u)
            for (w, uu), idx_u in edge_to_pipe_idx.items():
                if uu == u and pipe_enabled[idx_u] and codes[idx_u] > codes[idx_d]:
                    codes[idx_d] = codes[idx_u]
                    changed = True
        if not changed:
            break
    for i, c in enumerate(codes):
        if pipe_enabled[i] and c >= 0:
            pipe_diameters[i] = weight_list[c]
            
def quantize_lid_vars_4levels(lid_vars: np.ndarray) -> np.ndarray:
    """对 shape=(2*S,) 的 [BC1,PP1, BC2,PP2,...] 先量化到四档，再修复 BC+PP ≤ 0.10。"""
    def _quantize_to_4(v: float) -> float:
        v = float(np.clip(v, 0.0, 0.10))
        if v < 0.0125:      return 0.0
        elif v < 0.0375:    return 0.025
        elif v < 0.075:     return 0.05
        else:               return 0.10

    
    lv = lid_vars.astype(float).copy()
    for i in range(0, len(lv), 2):
        bc = _quantize_to_4(lv[i])
        pp = _quantize_to_4(lv[i+1])
        if bc + pp > 0.10 + 1e-12:
            # 保留较大者，把较小者降到不超过剩余额度的最大档
            if bc >= pp:
                remain = max(0.0, 0.10 - bc)
                pp = max([a for a in LID_LEVELS if a <= remain] or [0.0])
            else:
                remain = max(0.0, 0.10 - pp)
                bc = max([a for a in LID_LEVELS if a <= remain] or [0.0])
        lv[i], lv[i+1] = bc, pp
    return lv
            
def vector_to_design(
    x: np.ndarray,
    N_pipes: int,
    M_nodes: int,
    weight_list: List[float],
    *,
    enforce_monotone: bool,
    edge_to_pipe_idx: Optional[Dict[Tuple[int,int], int]] = None,
    fixed_outlet_mask: Optional[np.ndarray] = None,  # 若要和 GA 一样固定系统末端，可传入形如( M_nodes,)的0/1掩码
    #mst_indices: Optional[Sequence[int]] = None,
)-> tuple[np.ndarray, list[Optional[float]], np.ndarray, np.ndarray]:
    """
    纯解码器：把 BO 连续向量 x → (pipe_enabled, pipe_diameters, outlet_bits, lid_vars)
    约定：x = [states_norm( N_pipes ), outlet_bits( M_nodes ), lid_vars( 2*S )]
         其中 states_norm ∈ [0,1]，表示“直径档位/L”，不是 0/1！
    """
    L = len(weight_list)
    assert x.shape[0] >= N_pipes + M_nodes, "x 维度不足"

    # 1) 直径档位解码（注意：不再使用 *10 的硬编码）
    states_norm = np.clip(x[:N_pipes], 0.0, 1.0)
    states = np.clip(np.round(states_norm * L), 0, L).astype(int)  # 0..L
    pipe_enabled = (states > 0).astype(int)

    pipe_diameters: List[Optional[float]] = [None] * N_pipes
    for i, s in enumerate(states):
        if s >= 1:
            pipe_diameters[i] = weight_list[s - 1]

    # 2) 出口位
    outlet_slice = x[N_pipes : N_pipes + M_nodes]
    outlet_bits = np.round(np.clip(outlet_slice, 0.0, 1.0)).astype(int)
    if fixed_outlet_mask is not None:
        # 与 GA 对齐：若要求系统末端必须为 1，则覆盖/或取最大
        outlet_bits = np.maximum(outlet_bits, fixed_outlet_mask.astype(int))
        # 如果严格只允许固定的为 1，可改为：outlet_bits = fixed_outlet_mask.astype(int).copy()

    # 3) LID 变量（BC,PP, BC,PP, ...）
    lid_vars = x[N_pipes + M_nodes : ].astype(float).copy()
    lid_vars = quantize_lid_vars_4levels(lid_vars)

    # 4) （可选）全网单调修复：上游 ≤ 下游
    if enforce_monotone and edge_to_pipe_idx is not None:
        # 需要你保留之前给过的 repair_monotone_global(...)
        repair_monotone_global(pipe_enabled, pipe_diameters, weight_list, edge_to_pipe_idx, n_iter=10)

    return pipe_enabled, pipe_diameters, outlet_bits, lid_vars

'''
def _sample_lid_pair_continuous() -> Tuple[float, float]:
    """连续 LID：BC~U[0,0.1]，PP~U[0,0.1-BC]，保证两者之和≤0.1"""
    bc = np.random.rand() * 0.10
    pp = np.random.rand() * (0.10 - bc)
    return bc, pp
'''


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

'''
def lhs_initial_population(N_init, D, max_try=100):
    """
    用 Latin Hypercube Sampling (LHS) 初始化结构可行解向量。
    """
    raw_lhs = lhs(n=D, samples=N_init, criterion='maximin')
    samples = []
    for i in range(N_init):
        for attempt in range(max_try):
            try:
                # 将 LHS 采样点通过已有结构初始化函数转换为结构向量
                x_decoded = init_population_feasible_bo(
                    N_pipes=N_pipes,
                    M_nodes=M,
                    mst_indices=sorted(mst_indices),
                    edge_to_pipe_idx=edge_to_pipe_idx,
                    node_to_idx=node_to_idx,
                    system_outfall_nodes=system_outfall_node,
                    weight_list=weight_list,
                    catchment_areas=catchment_areas,
                    add_prob=0.15,
                    rng=None,
                    max_try=200
                )[-1]  # 只取设计向量 x

                if x_decoded is not None and isinstance(x_decoded, np.ndarray) and x_decoded.ndim == 1:
                    samples.append(x_decoded)
                    break
            except Exception as e:
                print(f"Attempt {attempt} failed for sample {i}: {e}")
                continue
    return samples
'''

def init_population_feasible_bo(
    N_pipes: int,
    M_nodes: int,
    mst_indices: List[int],                             # 骨架（MST）边对应的 pipe 索引
    edge_to_pipe_idx: Dict[Tuple[int, int], int],       # (u,v)->idx
    node_to_idx: Dict[int, int],                        # 节点ID->0..M-1（用于出口位）
    system_outfall_nodes: List[int],                    # 系统末端节点列表（固定为出口）
    weight_list: List[float],                           # 可选直径档（按小到大）
    catchment_areas: List[float],                       # 子汇面积（长度=子汇数S）
    add_prob: float = 0.15,                             # 非骨架随机加边概率（对齐 GA）
    rng: Optional[np.random.Generator] = None,
    max_try: int = 200
):
    """
    一步生成‘可行个体’，等价 GA 的 init_population_feasible：
      - 启用边：MST 必开 + 非骨架以 add_prob 打开
      - 管径：对启用边随机给档位，再做全网单调修复（上游≤下游）
      - 出口：固定系统末端为 1，其它为 0
      - LID：按 lid_mode 生成 (BC,PP)，并确保每个子汇 BC+PP≤0.10
    返回：
      pipe_enabled: (N_pipes,) int {0,1}
      pipe_diameters: List[float or None]
      outlet_bits: (M_nodes,) int {0,1}
      lid_vars: (2*S,) float  [BC1,PP1, BC2,PP2, ...]
      states_norm: (N_pipes,) float  # 档位归一化到[0,1]供 BO 向量化
      x: (N_pipes + M_nodes + 2*S,)  # 设计向量（和 vector_to_design 约定一致）
    """
    if rng is None:
        rng = np.random.default_rng()
    L = len(weight_list)

    for _ in range(max_try):
        
        # 1) 启用边：MST 必开 + 非骨架以 add_prob 打开
        pipe_enabled = np.zeros(N_pipes, dtype=int)
        if len(mst_indices) > 0:
            pipe_enabled[np.array(mst_indices, dtype=int)] = 1
        flip = rng.random(N_pipes) < add_prob
        pipe_enabled = np.where((pipe_enabled == 1) | flip, 1, 0).astype(int)

        # 2) 初始管径档位：对启用边随机在 [1..L] 采样，禁用边为 0 档（占位）
        states = np.zeros(N_pipes, dtype=int)
        k = np.where(pipe_enabled == 1)[0]
        if k.size > 0:
            states[k] = rng.integers(low=1, high=L + 1, size=k.size, endpoint=False) + 0  # [1..L]
        # 归一化到 [0,1] 给 BO 用；注意后面会再修复，但归一化的向量也要同步更新
        states_norm = states / float(L)

        # 3) 出口位：固定掩码或由节点列表生成
        if OUTLET_FIXED_MASK is not None:
            outlet_bits = OUTLET_FIXED_MASK.astype(int).copy()
            assert outlet_bits.shape == (M_nodes,)
        else:
            outlet_bits = np.zeros(M_nodes, dtype=int)
            for n in system_outfall_nodes:
                j = node_to_idx.get(n, None)
                if j is not None: outlet_bits[j] = 1

        # 4) LID：连续或 GA 离散
        lid_pairs = []
        '''
        sample_pair = _sample_lid_pair_continuous 
        for _ in range(len(catchment_areas)):
            lid_pairs.append(sample_pair())
        # 展平为 [BC1, PP1, BC2, PP2, ...]
        lid_vars = np.array([v for pair in lid_pairs for v in pair], dtype=float)
        '''
        for _ in range(len(catchment_areas)):
            bc = rng.choice(LID_LEVELS)
            pp = rng.choice(LID_LEVELS)
            if bc + pp > 0.10 + 1e-12:
                if bc >= pp:
                    remain = max(0.0, 0.10 - bc)
                    pp = max([a for a in LID_LEVELS if a <= remain] or [0.0])
                else:
                    remain = max(0.0, 0.10 - pp)
                    bc = max([a for a in LID_LEVELS if a <= remain] or [0.0])
            lid_pairs.append((bc, pp))
        lid_vars = np.array([v for pair in lid_pairs for v in pair], dtype=float)
        
        # 5) 把 states → 直径，并做全网单调修复
        pipe_diameters: List[Optional[float]] = [None] * N_pipes
        for i, s in enumerate(states):
            if s >= 1:
                pipe_diameters[i] = weight_list[s - 1]
        repair_monotone_global(pipe_enabled, pipe_diameters, weight_list, edge_to_pipe_idx, n_iter=10)

        # 6) 生成 BO 用设计向量 x（注意：直径档位部分是‘档位比例’而不是 0/1！）
        x = np.concatenate([states_norm.astype(np.float32),
                            outlet_bits.astype(np.float32),
                            lid_vars.astype(np.float32)], axis=0)
        
        # 结构约束：兼容 bool 或 (bool, msg)
        pipe_enabled, pipe_diams, outlet_bits, lid_vars = vector_to_design(x,
                                                                            N_pipes=N_pipes,
                                                                            M_nodes=M,
                                                                            weight_list=weight_list,
                                                                            enforce_monotone=True,
                                                                            edge_to_pipe_idx=edge_to_pipe_idx,
                                                                            fixed_outlet_mask=OUTLET_FIXED_MASK,  # ★ 关键
                                                                            #mst_indices=mst_indices
                                                                            )
        idx_to_node = {v: k for k, v in node_to_idx.items()}
        outlet_idx = np.flatnonzero(np.asarray(outlet_bits)).tolist()
        outlet_nodes_list = [idx_to_node[j] for j in outlet_idx]  # ← 这是 Python 列表

        _struct_ret = check_structure_validity(pipe_enabled, pipe_diams, outlet_nodes_list,
                                               edge_to_pipe_idx, G_full,
                                               subcatch_outlets=None,
                                               require_single_component=True,     # 整网一个弱连通组件
                                               require_reachability=False          # 每个节点都能到至少一个出口
                                               )
        if isinstance(_struct_ret, tuple):
            ok_struct = bool(_struct_ret[0])
            msg_s = _struct_ret[1] if len(_struct_ret) > 1 else ""
        else:
            ok_struct = bool(_struct_ret)
            msg_s = ""

        # LID 约束：兼容不同签名 & 返回形式
        try:
            _lid_ret = check_lid_constraints(lid_vars, catchment_areas, max_frac=0.1)
        except TypeError:
            try:
                _lid_ret = check_lid_constraints(lid_vars, catchment_areas)
            except TypeError:
                _lid_ret = check_lid_constraints(lid_vars)

        if isinstance(_lid_ret, tuple):
            ok_lid = bool(_lid_ret[0])
            msg_l = _lid_ret[1] if len(_lid_ret) > 1 else ""
        else:
            ok_lid = bool(_lid_ret)
            msg_l = ""
        if ok_struct and ok_lid:
            return pipe_enabled, pipe_diams, outlet_bits, lid_vars, states_norm, x
        
        fail_stats = defaultdict(int)
        if not ok_struct or not ok_lid:
            tag_s = "struct_ok" if ok_struct else "struct_ng"
            tag_l = "lid_ok"    if ok_lid    else "lid_ng"
            fail_stats[(tag_s, tag_l)] += 1
            continue
        print("[init_population_feasible_bo] fail stats:", dict(fail_stats))
    raise RuntimeError("init_population_feasible_bo: cannot find feasible sample within max_try")

def generate_inp_file(template_inp_path, output_inp_path,
                      pipe_enabled, pipe_diameters, conduits_info,
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
        if to in system_outfall_node:
            if to in used_outfalls:
                # 禁用掉第二条以上接到同一节点的管线
                continue
            # 第一次遇到，记录
            used_outfalls.add(to)

        # 3) 同时，确保这个 outfall 也不会被当作 from_node
        if frm in system_outfall_node:
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
        
        '''
        # if it is executable, comment it out
        if d is None or d <= 0:
            # print when d is None or <= 0
            print(f"Warning: Pipe {pipe_id} has invalid diameter {d}. Using default 0.25 m.")
            d = 0.25
        '''
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
        if node in  system_outfall_node:
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

def project_to_feasible_continuous(
    x_raw: np.ndarray,
    *,
    N_pipes: int, M_nodes: int, L: int,
    node_list, conduits_info,
    mst_indices, system_outfall_nodes,
    catchment_areas, max_lid_frac=0.10,
    min_level=1, force_level='mid'  # 'min'|'mid'|'max'
):
    """
    连续编码版的可行域投影：
      - 把必须启用的边（MST、连通所需）把 states_norm 提升到 >= min_level/L；
      - 如需更稳，可把这些边直接设为 'max'（=1.0），与老版“1→10档”对应；
      - 出口从 system_outfall_nodes 里保证 2 个；
      - LID 做上限缩放或 4 档量化+修复。
    返回 x_fixed。
    """
    x = np.array(x_raw, dtype=float).copy()
    assert x.ndim == 1

    # 1) 切段
    states_norm = x[:N_pipes]                    # 连续档位比例
    out_bits    = x[N_pipes:N_pipes+M_nodes]
    lid_vars    = x[N_pipes+M_nodes:]

    # 2) 出口：保证是 outfall 节点中的两个
    chosen = [node_list[i] for i,b in enumerate(out_bits>0.5) if b]
    # 补足/裁剪
    chosen = list(dict.fromkeys([n for n in chosen if n in system_outfall_nodes]))[:2]
    import random
    while len(chosen) < 2:
        c = random.choice(system_outfall_nodes)
        if c not in chosen:
            chosen.append(c)
    out_bits_fixed = np.array([1 if n in chosen else 0 for n in node_list], dtype=float)

    # 3) 必须启用的边集合：MST + 保证各出口与 chosen[0] 连通的路径
    import networkx as nx
    G = nx.DiGraph()
    edge_to_idx = {}
    for i, row in enumerate(conduits_info):
        u, v = int(row[1]), int(row[2])
        G.add_edge(u, v)
        edge_to_idx[(u, v)] = i

    must_on = set(mst_indices)
    base = chosen[0]
    for o in chosen:
        if nx.has_path(G, o, base):
            path = nx.shortest_path(G, o, base)
            for u, v in zip(path, path[1:]):
                idx = edge_to_idx[(u, v)]
                must_on.add(idx)

    # 4) 把这些“必须开启的边”的 states_norm 提升到 >0
    if force_level == 'max':
        target = 1.0
    elif force_level == 'mid':
        target = (min_level + L) / (2*L)  # 粗略中位比例
    else:  # 'min'
        target = max(min_level / L, 1e-6)

    states_norm_fixed = states_norm.copy()
    for idx in must_on:
        states_norm_fixed[idx] = max(states_norm_fixed[idx], target)

    # 5) LID：按 cap 缩放（或四档量化再修复 ≤0.10）
    # 5.1 直接缩放版：
    lid = lid_vars.copy()
    for i in range(len(catchment_areas)):
        bc, pp = lid[2*i], lid[2*i+1]
        tot = bc + pp
        if tot > max_lid_frac:
            fac = max_lid_frac / tot
            lid[2*i] *= fac; lid[2*i+1] *= fac

    # 如果你已经统一成四档，可以替换为：
    # lid = quantize_lid_vars_4levels(lid)  # 并包含 ≤0.10 修复

    # 6) 写回
    x_fixed = x.copy()
    x_fixed[:N_pipes] = np.clip(states_norm_fixed, 0.0, 1.0)
    x_fixed[N_pipes:N_pipes+M_nodes] = out_bits_fixed
    x_fixed[N_pipes+M_nodes:] = lid
    return x_fixed


twhole = time.time()  # [TIMING] total time    
class MixedKernel(gpytorch.kernels.Kernel): # adjust variables********/figure out the meaning of kernel
    def __init__(self, cont_kernel, disc_kernel):
        super(MixedKernel, self).__init__()
        self.cont_kernel = cont_kernel
        self.disc_kernel = disc_kernel
    
    def forward(self, x1, x2, diag=False, **params):
        x1_cont = x1[..., :d_cont]
        x2_cont = x2[..., :d_cont]
        x1_disc = x1[..., d_cont:]
        x2_disc = x2[..., d_cont:]
        k_cont = self.cont_kernel(x1_cont, x2_cont, diag=diag, **params)
        k_disc = self.disc_kernel(x1_disc, x2_disc, diag=diag, **params)
        return k_cont * k_disc

# Sub-kernels:
cont_kernel = gpytorch.kernels.RBFKernel(ard_num_dims=d_cont)
disc_kernel = gpytorch.kernels.RBFKernel(ard_num_dims=d_disc)
combined_kernel = gpytorch.kernels.ScaleKernel(MixedKernel(cont_kernel, disc_kernel))

class MixedGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, input_dim):
        super(MixedGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = combined_kernel
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
    
##########################
# Reorder the Design Vectors for GP Input
##########################
def reorder_for_kernel(x, N_pipes: int = N_pipes, M_nodes: int = M, S: int = 18):
    X = np.asarray(x, dtype=np.float32)
    single = (X.ndim == 1)
    if single: X = X[None, :]

    # 连续：前 N_pipes
    states_norm = X[:, :N_pipes]

    # 离散：出口位二值化
    outlet_bits = (X[:, N_pipes:N_pipes+M_nodes] > 0.5).astype(np.float32)

    # 离散：LID 四档量化（调用统一函数）
    lid = X[:, N_pipes+M_nodes : N_pipes+M_nodes + 2*S]
    lid_q = np.vstack([quantize_lid_vars_4levels(row) for row in lid])

    X_proj = np.concatenate([states_norm, outlet_bits, lid_q], axis=1)
    return X_proj[0] if single else X_proj


def reorder_batch(X):
    return np.array([reorder_for_kernel(x) for x in X])

# Define the Expected Improvement (EI) acquisition function.

def expected_improvement(X, f_best, model, xi=0.0):
    """
    X: (n, D) torch.float，已经 reorder 过
    f_best: float，当前最小的目标值（极小化）
    model: 训练好的 GP
    xi: 探索系数（0~0.05 常用；越大越探索）
    """
    model.eval()
    with torch.no_grad():
        post = model(X)
        mu   = post.mean            # (n,)
        var  = post.variance.clamp_min(1e-12)
        sigma = var.sqrt()

        # 极小化：imp = f_best - mu - xi
        imp = f_best - mu - xi
        Z   = imp / sigma

        # EI = imp * Phi(Z) + sigma * phi(Z)；sigma≈0 时令 EI=0
        normal = torch.distributions.Normal(0.0, 1.0)
        ei  = imp * normal.cdf(Z) + sigma * torch.exp(-0.5 * Z**2) / np.sqrt(2*np.pi)
        ei = torch.where(sigma <= 1e-12, torch.zeros_like(ei), ei)
    return ei

# if I want to use a piecewise schedule for xi in BO, if no need, just commend out
def xi_schedule_piecewise(it, N_initial, max_iter,
                          xi_hi=0.1, xi_mid=0.05, xi_lo=0.01,
                          p1=0.33, p2=0.67):
    """it 是当前迭代编号（从 N_initial 开始）；返回本轮 xi。"""
    k = max(0, it - N_initial)
    p = k / max(1, max_iter)  # 进度 0~1
    if p < p1:
        return xi_hi
    elif p < p2:
        return xi_mid
    else:
        return xi_lo

def Bayesian_Opt(inp_file, rpt_file):
    samples = []
    q_init = []
    c_init = []
    #rpt_files = []
    tgenersample = time.time()  # [TIMING] sampling
    
    for i in range(N_initial):
        X_init = init_population_feasible_bo(
            N_pipes=N_pipes,
            M_nodes=M,
            mst_indices=sorted(mst_indices),
            edge_to_pipe_idx=edge_to_pipe_idx,
            node_to_idx=node_to_idx,
            system_outfall_nodes=system_outfall_node,
            weight_list=weight_list,
            catchment_areas=catchment_areas,
            add_prob=0.15,
            rng=None,
            max_try=200
        )[-1]  # 只取设计向量 x
        if X_init is None or not isinstance(X_init, np.ndarray) or X_init.ndim != 1:
            print(f"Sample {i} returned bad x (None or not 1D), skipping:", X_init)
            continue
        if np.any(np.isnan(X_init)) or np.any(np.isinf(X_init)):
            print(f"Sample {i} contains NaN/Inf, skipping.")
            continue
        X_init = project_to_feasible_continuous(
            X_init,
            N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
            node_list=node_list, conduits_info=conduits_info,
            mst_indices=sorted(mst_indices),
            system_outfall_nodes=system_outfall_node,
            catchment_areas=catchment_areas,
            max_lid_frac=0.10,
            force_level='mid'     # ← 与你候选/评估处保持一致
        )
        
        pipe_enabled, pipe_diams, outlet_bits, lid_vars = vector_to_design(#X_init[0],   # 只取第一个样本
                                                                            X_init,
                                                                            N_pipes=N_pipes,
                                                                            M_nodes=M,
                                                                            weight_list=weight_list,
                                                                            enforce_monotone=True,
                                                                            edge_to_pipe_idx=edge_to_pipe_idx,
                                                                            fixed_outlet_mask=OUTLET_FIXED_MASK,  # ★ 关键
                                                                            #mst_indices=mst_indices
                                                                            )

        generate_inp_file(template_inp_path, inp_file, pipe_enabled, pipe_diams, conduits_info, lid_vars, loops, subcatch_info, coords_info)    
        samples.append(X_init)
    
    #samples = lhs_initial_population(N_init=N_initial, D=D)
        # check if all samples are the same, if so, print a warning
        if len(samples) > 1 and np.allclose(samples[-1], samples[0]):
            print(f"Warning: Sample {i} is the same as the first sample, all samples so far are identical.")
            #break  # Uncomment to stop if all samples are the same
        if len(samples) == 0:
            raise RuntimeError("Failed to generate a valid initial design.")    
        X_init = np.stack(samples, axis=0).astype(np.float32)
        print(f"Generated {len(samples)} valid samples in {time.time() - tgenersample:.2f} seconds.")
        np.random.shuffle(X_init)
        idx = np.random.permutation(len(X_init))
        tsampcheck = time.time()                         # [TIMING] sample check
    
        #conduits_mod = conduits_info
        #coords_mod   = coords_info
        
        #print(f">>> Running sample {i} …", flush=True)

        #rpt_files.append(rpt_file)
        c = compute_LCC(inp_file,cost_per_diameter,lid_vars,
                            cost_lid_bc=200, cost_lid_pp=150,
                            om_lid_bc=0.08, om_lid_pp=0.04,
                            discount_rate=0.02, project_lifetime=30)
        print(f"    <- SWMM subprocess finished")
        # Define objective as a combination. For example: y = qmax + alpha * lcc, alpha=1 for simplicity.
        # qmax 
        subprocess.run(
            [swmm_exe_path, inp_file, rpt_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        q = get_system_qmax_from_rpt(rpt_file)
        
        q_init.append(q)
        c_init.append(c)
        
        
    q_arr = np.array(q_init, dtype=np.float32)
    c_arr = np.array(c_init, dtype=np.float32)
    y_init = []
    qmax_ini_list = []
    lcc_ini_list = []
    for q, c in zip(q_arr, c_arr):
        qn_ini = (q - Q_MAX_MIN) / (Q_MAX_MAX - Q_MAX_MIN + EPS)
        cn_ini = (c - LCC_MIN) / (LCC_MAX - LCC_MIN + EPS)
        #qn_ini = float(np.clip(qn_ini, 0.0, 1.0))
        #cn_ini = float(np.clip(cn_ini, 0.0, 1.0))
        y_init.append(w2 * qn_ini + w1 * cn_ini)
        qmax_ini_list.append(qn_ini)
        lcc_ini_list.append(cn_ini)
        # print lcc and qmax, separate them
        print(f"Sample: qmax={q:.2f}, lcc={c:.2f}, normalized y={y_init[-1]:.4f}, Q_MAX_MIN={Q_MAX_MIN:.4f}, Q_MAX_MAX={Q_MAX_MAX:.4f}, LCC_MIN={LCC_MIN:.4f}, LCC_MAX={LCC_MAX:.4f}")

    X_init = X_init.astype(np.float32)
    X_init = X_init[idx]
    y_init = np.array(y_init, dtype=np.float32)
    y_init = y_init[idx]
    qmax_ini_list = np.array(qmax_ini_list, dtype=np.float32)
    lcc_ini_list = np.array(lcc_ini_list, dtype=np.float32)
    qmax_ini_list = qmax_ini_list[idx]
    lcc_ini_list  = lcc_ini_list[idx]
    print(f"[init {i}] total: {time.time()-tsampcheck:.3f}s")       # [TIMING]
    print("=== Debug: y_init summary ===")
    print("  y_init shape:", y_init.shape)
    print("  unique values in y_init:", np.unique(y_init))
    print("  min, max, mean:", y_init.min(), y_init.max(), y_init.mean())
    
    # ========================
    # Bayesian Optimization Loop
    # ========================
    times   = [] # for the time plot
    t_train_list = []
    t_cand_list  = []
    t_ei_list    = []
    t_feas_list  = []
    t_inp_list   = []
    t_swmm_list  = []
    t_qmax_list = []
    t_obj_list   = []
    t_iter_list  = []

    lcc_list = []
    qmax_list = []
    y_next_list = []
    #y_bo = []
    iter_time = time.time()
    
    for it in range(N_initial, N_initial + max_iter):
        t_iter_start = time.time()                                        # [TIMING] 单次迭代总耗时
        # a) Reorder and convert to torch for GP
        t_reorder_start = time.time()
        #X_init = reorder_batch(X_init)
        t_reorder = time.time() - t_reorder_start
        
        t_train_start = time.time() # [TIMING] GP 训练耗时
        print(">>> Training GP surrogate model")
        # Convert initial data to torch tensors.
        #train_x = torch.tensor(X_init, dtype=torch.float)
        train_x = torch.tensor(reorder_batch(X_init), dtype=torch.float)  # reorder for GP
        train_y = torch.tensor(y_init, dtype=torch.float)


        #b) Train the GP Surrogate with Combined Kernel
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        model = MixedGPModel(train_x, train_y, likelihood, input_dim=D)

        # Find optimal model hyperparameters by training the GP.
        model.train(); likelihood.train()
        optimizer_gp = torch.optim.Adam(model.parameters(), lr=0.001)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

        training_iter = 100 # 100 iterations for training the GP
        # the purpose of this is to optimize the hyperparameters of the GP model, minimize the loss
        for i in range(training_iter):
            optimizer_gp.zero_grad()
            output = model(train_x)
            loss = -mll(output, train_y)
            loss.backward()
            optimizer_gp.step()
            # Uncomment the next line to print loss occasionally
            # if (i+1) % 20 == 0:
            #     print(f"Iteration {i+1}/{training_iter} - Loss: {loss.item():.3f}")
            if (i+1) % 20 == 0 or i == 0:
                print(f"GP iter {i+1}/{training_iter} - Loss: {loss.item():.5f}")

        model.eval(); likelihood.eval()
        t_train = time.time() - t_train_start
        t_train_list.append(t_train)
        print(f"[Iter {it}] GP training: {t_train:.3f}s")
        
        # c) Propose next point via EI
        t_cand_start = time.time()                                            # [TIMING] EI 计算耗时
        f_best = train_y.min().item() # 取目前已评估样本 train_y（目标值）的最小值，作为当前最好值 𝑓∗（你是在做最小化）。.item() 把 0-维张量变成 Python 标量
        xi = xi_schedule_piecewise(it, N_initial, max_iter)
        # build hybrid candidate pool: random + local + top-hist
        n_random = 200
        n_local  = 500
        n_best   = 200
        best_x = X_init[np.argmin(y_init)]
        X_candidates = []
        for _ in range(n_random):
            x = np.random.rand(D)
            x = project_to_feasible_continuous(x, N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
                                               node_list=node_list, conduits_info=conduits_info,
                                               mst_indices=sorted(mst_indices),
                                               system_outfall_nodes=system_outfall_node,
                                               catchment_areas=catchment_areas,
                                               max_lid_frac=0.10, force_level='mid')
            X_candidates.append(x)
        for _ in range(n_local):
            noise = np.random.normal(loc=0.0, scale=0.05, size=D)
            x = np.clip(best_x + noise, 0.0, 1.0)
            x = project_to_feasible_continuous(x, N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
                                               node_list=node_list, conduits_info=conduits_info,
                                               mst_indices=sorted(mst_indices),
                                               system_outfall_nodes=system_outfall_node,
                                               catchment_areas=catchment_areas,
                                               max_lid_frac=0.10, force_level='mid')
            X_candidates.append(x)
        X_candidates.extend(X_init[np.argsort(y_init)[:n_best]])
        X_candidates = np.unique(np.round(X_candidates, 6), axis=0)
        
        # === Compute EI [MODIFIED]
        X_cand_torch = torch.tensor(reorder_batch(X_candidates), dtype=torch.float)
        ei_vals = expected_improvement(X_cand_torch, f_best, model, xi=xi)

        # Soft penalty for proximity to last_x
        last_x = reorder_for_kernel(X_init[-1])
        ei_np = ei_vals.detach().cpu().numpy().ravel()
        distances = np.linalg.norm(reorder_batch(X_candidates) - last_x, axis=1)
        ei_np -= 1.0 * np.exp(-distances**2 / 1e-5)  # soft repulsion
        x_next_raw = X_candidates[np.argmax(ei_np)]
        
        # === [PATCH] 轻量“邻域拣选”：仅用 GP 选择一个邻居（仍只跑 1 次 SWMM） ===
        k_local = 64
        neigh = []
        for _ in range(k_local):
            z = np.clip(x_next_raw + np.random.normal(0, 0.02, size=x_next_raw.shape), 0, 1)
            z = project_to_feasible_continuous(
                z, N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
                node_list=node_list, conduits_info=conduits_info,
                mst_indices=sorted(mst_indices),
                system_outfall_nodes=system_outfall_node,
                catchment_areas=catchment_areas,
                max_lid_frac=0.10, force_level='mid'
            )
            neigh.append(z)

        # 去重，避免全被投影到同一点
        neigh = np.unique(np.round(np.asarray(neigh), 6), axis=0)
        if neigh.shape[0] == 0:
            x_eval = x_next_raw.copy()
        else:
            with torch.no_grad():
                Nt = torch.tensor(reorder_batch(neigh), dtype=torch.float)
                post = model(Nt)
                mu = post.mean
                sigma = post.variance.clamp_min(1e-12).sqrt()
                beta = 1.5
                # 最小化：选 mu - beta * sigma 最小者（UCB for minimization）
                ucb = mu - beta * sigma
                x_eval = neigh[int(torch.argmin(ucb))]
        
        # d) Evaluate x_next via SWMM
        # 修复到可行
        x_next = project_to_feasible_continuous(x_eval,
                                                N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
                                                node_list=node_list, conduits_info=conduits_info,
                                                mst_indices=mst_indices,
                                                system_outfall_nodes=system_outfall_node,
                                                catchment_areas=catchment_areas,
                                                max_lid_frac=0.10,
                                                force_level='mid'   # 若要完全复刻老版“启用=最大档”的行为
                                                )
        pipe_enabled, pipe_diameters, outlet_selection, lid_vars = vector_to_design(x_next,
                                                                                    N_pipes=N_pipes,
                                                                                    M_nodes=M,
                                                                                    weight_list=weight_list,
                                                                                    enforce_monotone=True,
                                                                                    edge_to_pipe_idx=edge_to_pipe_idx,
                                                                                    fixed_outlet_mask=OUTLET_FIXED_MASK,
                                                                                    #mst_indices=mst_indices
                                                                                    )
        for idx in mst_indices:
            pipe_enabled[idx] = 1
            
        '''
        # 位→索引（或索引→节点ID）
        out_idx = np.flatnonzero(outlet_bits).tolist()
        # 如果 check_structure_validity 需要节点ID：
        idx_to_node = {v:k for k,v in node_to_idx.items()}
        outlet_nodes = [idx_to_node[i] for i in out_idx]
        
        if len(outlet_selection) == 0 or (not check_lid_constraints(lid_vars, catchment_areas)) or not (check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes,
                                                                                                                                 edge_to_pipe_idx=edge_to_pipe_idx, G_full=G_full, subcatch_outlets=None,
                                                                                                                                 require_single_component=True,     # 整网一个弱连通组件
                                                                                                                                 require_reachability=False)):
            y_next = 1e6
        else:
        '''
        t_inp_start = time.time()   # [TIMING] 生成输入文件耗时
        generate_inp_file(template_inp_path, inp_file, pipe_enabled, pipe_diameters,
                          conduits_info, lid_vars, loops, subcatch_info, coords_info)
        t_inp = time.time() - t_inp_start
        t_inp_list.append(t_inp)
        print(f"[Iter {it}] generate_inp_file: {t_inp:.3f}s")
        t_swmm_start = time.time()                                        # [TIMING] SWMM 运行耗时
        lcc = compute_LCC(inp_file,cost_per_diameter,lid_vars,
                            cost_lid_bc=200, cost_lid_pp=150,
                            om_lid_bc=0.08, om_lid_pp=0.04,
                            discount_rate=0.02, project_lifetime=30)
        #print(cost_per_diameter)
        #print(lcc)
        #lcc_list.append(lcc)
        t_swmm = time.time() - t_swmm_start
        t_swmm_list.append(t_swmm)
        print(f"[Iter {it}] SWMM run: {t_swmm:.3f}s")     # [TIMING]
        tqmax_start = time.time() 
        # qmax 
        subprocess.run(
            [swmm_exe_path, inp_file, rpt_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )# [TIMING] obtain Qmax
        qmax = get_system_qmax_from_rpt(rpt_file)
        #qmax_list.append(qmax)
        #print("print all qmax_list:", qmax_list)
        t_qmax = time.time() - tqmax_start
        t_qmax_list.append(t_qmax)
        print(f"[Iter {it}] parse RPT: {t_qmax:.3f}s")     # [TIMING]
        t_obj_start = time.time()
        qn = (qmax - Q_MAX_MIN) / (Q_MAX_MAX - Q_MAX_MIN + EPS) # ****to do: fixed qmax and qmin, unify norm range in both algorithms
        cn = (lcc  - LCC_MIN) / (LCC_MAX - LCC_MIN + EPS) # ****to do: same, and print out the specific values, find out ga curve (qmax), further optimize
        #qn = float(np.clip(qn, 0.0, 1.0))
        #cn = float(np.clip(cn, 0.0, 1.0))
        y_next = w2 * qn + w1 * cn
        
        t_obj = time.time() - t_obj_start
        t_obj_list.append(t_obj)
        print(f"[Iter {it}] objective: {t_obj:.3f}s (y={y_next:.4f}, qn={qn:.4f}, cn={cn:.4f})")
            
        y_next_list.append(y_next)
        lcc_list.append(cn)
        qmax_list.append(qn)
        #y_bo.append(y_next)
        times.append(time.time() - iter_time)
        
        # e) Append new sample
        X_init = np.vstack([X_init, x_next[None,:]])
        #y_total = np.concatenate([y_init, np.array([y_next], dtype=np.float32)])
        #y_init = np.array(list(y_init) + list(y_next_list), dtype=np.float32)
        y_init = np.append(y_init, np.float32(y_next))
        #lcc_ini_list  = np.array(list(lcc_ini_list)  + list(lcc_list),  dtype=np.float32)
        #qmax_ini_list = np.array(list(qmax_ini_list) + list(qmax_list), dtype=np.float32)
        qmax_ini_list = np.append(qmax_ini_list, np.float32(qn))
        lcc_ini_list  = np.append(lcc_ini_list,  np.float32(cn))

        t_iter = time.time() - t_iter_start
        t_iter_list.append(t_iter)
        print(f"[Iter {it}] TOTAL iteration time: {t_iter:.3f}s")
        #print(f"Iteration {it+1}: y_next = {y_next:.3f}, EI = {best_ei_value:.3e}")

        best_idx = np.argmin(y_init)
        #best_idx = np.argmin(y_bo) 
        #print("Best y =", y_init[best_idx], "at x =", X_init[best_idx])

        # Convert to a Python list if it's a numpy array
        y_history = list(y_init)
        #y_history = list(y_bo)  # Use the Bayesian optimization history

        # Compute the best-so-far sequence
        best_so_far = [min(y_history[:i+1]) for i in range(len(y_history))]
        #qmax_best_so_far = [min(qmax_total[:i+1]) for i in range(len(qmax_total))]
        #lcc_best_so_far = [min(lcc_total[:i+1]) for i in range(len(lcc_total))]
        #best_so_far = [min(y_history[:i+1]) for i in range(len(y_next_list))]
        qmax_best_so_far = [min(qmax_ini_list[:i+1]) for i in range(len(qmax_list))]
        lcc_best_so_far = [min(lcc_ini_list[:i+1]) for i in range(len(lcc_list))]

    print("Decoded outlet selection:", outlet_selection)
    print("Objective value at x_next:", y_next)
    # print out twhole
    print(f"Total time for iteration {it}: {time.time() - twhole:.3f}s")
    
    plt.subplot(1, 3, 1)
    plt.plot(np.arange(1, len(y_init)+1)[:N_initial],y_init[:N_initial], marker='o', linestyle='--', label='Initial 50 samples (y)')
    plt.plot(range(1, max_iter+1), best_so_far[N_initial:N_initial+max_iter+1],marker='o',label='BO per-iteration y')
    plt.xlabel('Iteration')
    plt.ylabel('Best Objective Value So Far')
    plt.title('Bayesian Optimization Convergence Curve')
    plt.tight_layout()

    iterations = np.arange(1, len(qmax_list) + 1)
    #iterations = np.arange(1, len(y_bo)+1)
    # lcc 曲线
    plt.subplot(1, 3, 2)
    plt.plot(iterations, lcc_best_so_far, marker='o', linestyle='-')
    plt.xlabel("Iteration")
    plt.ylabel("LCC (Cost)")
    plt.title("LCC Over Iterations")
    plt.grid(True)
    plt.tight_layout()


    # qmax 曲线
    plt.subplot(1, 3, 3)
    plt.plot(iterations, qmax_best_so_far, marker='o', linestyle='-')
    plt.xlabel("Iteration")
    plt.ylabel("qmax (L/s)")
    plt.title("qmax Over Iterations")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


    # obj values vs time
    plt.figure(figsize=(6,4))
    plt.plot(times, best_so_far, marker='o', linestyle='-')
    plt.xlabel("Elapsed Time (s)")
    plt.ylabel("Best Objective So Far")
    plt.title("Convergence vs Wall‐Clock Time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
        
if __name__ == "__main__":
    template_inp_path = "./AhvazNull.inp"
    swmm_exe_path = "./swmm5.exe" 
    inp_file = "BS_opt.inp"
    rpt_file = "BS_opt.rpt"
    best_solution = Bayesian_Opt(inp_file, rpt_file)


'''
template_inp_path = "./AhvazNull.inp"
swmm_exe_path = "./swmm5.exe" 
inp_file = "BS_opt.inp"
rpt_file = "BS_opt.rpt"
samples = []
for i in range(N_initial):
    X_init = init_population_feasible_bo(
        N_pipes=N_pipes,
        M_nodes=M,
        mst_indices=sorted(mst_indices),
        edge_to_pipe_idx=edge_to_pipe_idx,
        node_to_idx=node_to_idx,
        system_outfall_nodes=system_outfall_node,
        weight_list=weight_list,
        catchment_areas=catchment_areas,
        add_prob=0.15,
        rng=None,
        max_try=200
    )[-1]  # 只取设计向量 x 
    if X_init is None or not isinstance(X_init, np.ndarray) or X_init.ndim != 1:
        print(f"Sample {i} returned bad x (None or not 1D), skipping:", X_init)
        continue
    if np.any(np.isnan(X_init)) or np.any(np.isinf(X_init)):
        print(f"Sample {i} contains NaN/Inf, skipping.")
        continue    
    samples.append(X_init)
    # check if all samples are the same, if so, print a warning
    if len(samples) > 1 and np.allclose(samples[-1], samples[0]):
        print(f"Warning: Sample {i} is the same as the first sample, all samples so far are identical.")
        #break  # Uncomment to stop if all samples are the same
if len(samples) == 0:
    raise RuntimeError("Failed to generate a valid initial design.")    
X_init = np.stack(samples, axis=0).astype(np.float32)
print(X_init)

pipe_enabled, pipe_diams, outlet_bits, lid_vars = vector_to_design(X_init[0],   # 只取第一个样本
                                                                                   N_pipes=N_pipes,
                                                                                   M_nodes=M,
                                                                                   weight_list=weight_list,
                                                                                   enforce_monotone=True,
                                                                                   edge_to_pipe_idx=edge_to_pipe_idx,
                                                                                   fixed_outlet_mask=OUTLET_FIXED_MASK   # ★ 关键
                                                                                )

generate_inp_file(template_inp_path, inp_file, pipe_enabled, pipe_diams,conduits_info, lid_vars,loops, subcatch_info, coords_info)
y = compute_LCC(inp_file, cost_per_diameter, lid_vars,
                  cost_lid_bc=200, cost_lid_pp=150,
                  om_lid_bc=0.08, om_lid_pp=0.04,
                  discount_rate=0.02, project_lifetime=30)
# qmax 
subprocess.run(
    [swmm_exe_path, inp_file, rpt_file],
    check=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
)
qmax = get_system_qmax_from_rpt(rpt_file)
print("Initial sample LCC:", y)
print("Initial sample Qmax:", qmax)
'''