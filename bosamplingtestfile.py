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
weight_to_state = {
    w: commercial_diameters.index(w) + 1
    for w in commercial_diameters
}

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
D = 65 + 18 + M + 36 # 设计向量总维度 (dynamic)
d_cont = 36
d_disc = 65 + 18 + M

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

# 取第一个作为主 outfall
#system_outfall_node = system_outfalls[0]
#print("Using system outfall node =", system_outfall_node)



def sample_weights_for_network(pipe_enabled, weights):
    """
    对启用的管道采样：——
      1）骨架（mst_pipe_order）上按拓扑保证 up≤down；
      2）非骨架边，在 [max_up, min_down] 范围内采样。
    返回长度 = 管道总数的列表，元素要么是 None（未启用），要么是直径值。
    """
    sampled = {}

    # 1) 在骨架 DAG 上按 mst_pipe_order 采样
    for p in mst_pipe_order:
        if not pipe_enabled[p]:
            continue
        else:
            sampled[p] = weights[0]
        # upstream preds & downstream succs on the MST
        ups   = skel_preds.get(p, [])  
        downs = [q for q in mst_pipe_order if p in skel_preds.get(q, [])]

        # lower bound = max weight among preds, or min-weight if none
        lower = max(sampled.get(u, weights[0]) for u in ups) if ups else weights[0]
        # upper bound = min weight among succs, or max-weight if none
        upper = min(sampled.get(d, weights[-1]) for d in downs) if downs else weights[-1]

        # pick a valid weight in [lower, upper]
        candidates = [w for w in weights if lower <= w <= upper]
        sampled[p] = random.choice(candidates) 

    # 2) 对非骨架但启用的管道也做上下界采样
    for i, en in enumerate(pipe_enabled):
        if not en or i in sampled:
            continue

        # 上界和下界从骨架上找
        # u→v 是这条管道的上下游节点
        u = int(conduits_info[i,1])
        v = int(conduits_info[i,2])

        # 找到所有骨架管道 j，其下游节点 == u → 取它们的 sampled[j] 作为下限候选
        ups = [j for j in mst_indices if int(conduits_info[j,2]) == u]
        lower = max(sampled[j] for j in ups) if ups else weights[0]

        # 找到所有骨架管道 k，其上游节点 == v → 取它们的 sampled[k] 作为上限候选
        downs = [k for k in mst_indices if int(conduits_info[k,1]) == v]
        upper = min(sampled[k] for k in downs) if downs else weights[-1]

        # 如果上下界反过来（upper < lower），就退化成下限
        if upper < lower:
            valid = [w for w in weights if w >= lower]
        else:
            valid = [w for w in weights if lower <= w <= upper]
        if not valid:
            # 最坏也只能用全局最小直径
            valid = [weights[0]]
        sampled[i] = random.choice(valid)

    # 3) 把 disabled 的管道设为 None
    return [ sampled[i] if pipe_enabled[i] else None
             for i in range(len(pipe_enabled)) ]
    
def vector_to_design(x):
    """
    从长度为 D = N_pipes + M + 36 的设计向量 x 中拆出：
      - pipe_enabled:  长度 N_pipes （0/1 启用）
      - pipe_diameters: 长度 N_pipes （连续值）
      - outlet_selection: 长度 M （0/1）
      - lid_vars: 长度 36 （连续值）
    """
    if not isinstance(x, np.ndarray) or x.ndim != 1:
        raise RuntimeError(f"[vector_to_design] expected 1-D array, got {x!r}")
    x = np.array(x, dtype=float)
    #N_pipes = conduits_info.shape[0]
    M       = len(node_list)     # 或者提前全局定义好


    # if x is nan, raise an error
    #if np.isnan(x).any():
    #    print(x)
    #    raise ValueError("Input vector x contains NaN values.")
    states = (x[0:N_pipes] * 10).round().astype(int)
    pipe_enabled   = (states > 0).astype(int)
    pipe_diameters = [
        (weight_list[s-1] if 1 <= s <= len(weight_list) else 0.0)
        for s in states
    ]

    # —— 出口选择位 —— 
    # 接着 N_pipes:N_pipes+M
    bits = x[N_pipes : N_pipes + M] > 0.5
    outlet_selection = [ node_list[i] for i,flag in enumerate(bits) if flag ]
    lid_vars = x[N_pipes + M : N_pipes + M + 36]

    return pipe_enabled, pipe_diameters, outlet_selection, lid_vars

        
## This is a draft of constraints (so far, grey only) implemented to BO
def check_structure_validity(pipe_enabled, pipe_diameters, outlet_nodes):
    """
    several constraints integrated together:
    1. number of outlet candidates >= 1
    2. all enabled pipes form a weakly connected graph (must be connected)
    3. upstream pipe diameter ≤ downstream pipe diameter
    """


    # constraint 1: number of outlet candidates >= 1
    #outlets = find_outlets(pipe_enabled)
    if not outlet_nodes:
        return False

    
    # constraint 2: all enabled pipes form a weakly connected graph 
    
    enabled_edges = [(u,v) for (u,v),idx in edge_to_pipe_idx.items() if pipe_enabled[idx]]
    G_enabled = G_full.edge_subgraph(enabled_edges).to_undirected()

    if not nx.is_connected(G_enabled):
        print("Warning: enabled pipes do not form a connected graph.")
        # raise error
        raise RuntimeError("Enabled pipes do not form a connected graph.")
        # return False
    else:
        print("Enabled pipes form a connected graph.")
    '''
    if not all_nodes_reach_outlet(pipe_enabled, outlet_nodes):
        return False
    '''
    
    # constraint 3: upstream pipe diameter ≤ downstream pipe diameter
    for (u,v), idx in edge_to_pipe_idx.items():
        if pipe_enabled[idx]:
            down = pipe_diameters[idx]
            # 遍历所有 u 的前驱 w→u
            for (w,uu), j in edge_to_pipe_idx.items():
                if uu == u and pipe_enabled[j]:
                    up = pipe_diameters[j]
                    if up and down and up > down:
                        return False

    return True  # if all constraints are satisfied


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
    
def sample_valid_design(conduits_info, catchment_areas, weight_list, weight_to_state, maxtry = 100):
    """
    Draw one x∈[0,1]^D that, after decoding, satisfies
       • len(outlet_selection)>0
       • check_lid_constraints(...)
       • check_structure_validity(...)
    """
    # add a counter to check which constraint is violated
    no_system_outfall = lid_fail = struct_fail = 0
    N = len(conduits_info) 
    for attempt in range(maxtry):
        # --- 步骤 1: 随机打开/关闭管道骨架外的边，骨架上的必须打开 ---
        pipe_enabled = [
            1 if i in mst_indices else random.choice([0,1])
            for i in range(N)
        ]
        
        
        # 2) 拓扑序+权重采样：保证每条启用管道的 downstream_weight >= 所有直接 upstream_weight
        sampled_weights = sample_weights_for_network(
            pipe_enabled,
            weight_list
        )
        # sampled_weights: dict {pipe_idx: weight_value}
        
        # 3) Pipe states (65 ints in [0,10])                                      
        pipe_states = []
        for en, w in zip(pipe_enabled, sampled_weights):
            if not en:
                pipe_states.append(0)
            else:
                pipe_states.append(weight_to_state[w])
        
        # 4) dynamic test: ensure pre‑picked system outfalls are still true ends
        # build subgraph of enabled edges
        enabled_edges = [(u,v) for (u,v), idx in edge_to_pipe_idx.items() if pipe_enabled[idx]]
        G_sub = nx.DiGraph(); G_sub.add_edges_from(enabled_edges)
        # find degree‑0 nodes not in subcatch (true network ends)
        ends = [n for n in G_sub.nodes()
                if G_sub.out_degree(n)==0 and n not in subcatch_outlets]
        # test validity
        if len(ends) < 2 or any(o not in ends for o in system_outfall_node):
            no_system_outfall += 1
            continue
        sel = system_outfall_node
        
        # 5) outlet_bits: mark exactly the two chosen system outfalls
        outlet_bits = [1 if n in system_outfall_node else 0 for n in node_list]
             
        
        # 7) LID ratios: for each of 18 subcatchments,
        
        #    sample BC∈[0,0.1], then PP∈[0, 0.1−BC]
        raw_bc = []
        raw_pp = []
        for A in catchment_areas:
            cap = 0.1 * A
            if cap <= 0:
                # no area ⇒ force zero LID
                raw_bc.append(0.0)
                raw_pp.append(0.0)
            else:
                b = random.uniform(0, cap)
                p = random.uniform(0, cap - b)
                raw_bc.append(b / A)
                raw_pp.append(p / A)

            
            
        # 1) 先动态取管道数和出口数
        #N_pipes = conduits_info.shape[0]
        # node_list, M 也需要是全局或传参的
        M = len(node_list)  

        # 2) 设计向量长度
        D = N_pipes + M + 36    
        
        # 3) 初始化
        x = np.zeros(D, dtype=np.float32)
        #print(x) 
        # pipe_states → 0–64
        try:
            x[0:N_pipes] = np.array(pipe_states)/10.0
            #print("piepe_states executed")
        except Exception as e:
            print(f"Error setting pipe states: {e}")
            continue
        
        x[N_pipes : N_pipes + M] = outlet_bits
        if np.isnan(x).any():
            print("DEBUG: x is nan when outlet", outlet_bits)
            continue
        #else:
        #    print("outlet_bits executed")
        # LID → x[65+M : 65+M+36]
        lid_flat = [val for pair in zip(raw_bc, raw_pp) for val in pair]
        start = N_pipes + M
        x[start : start + 36] = np.array(lid_flat, dtype=np.float32)
        if np.isnan(x).any():
            print("DEBUG: x is nan when LID", lid_flat)
            continue
        #else:
        #    print("LID executed")
        
        # print  x
        #print(f"DEBUG: Attempt {attempt+1}, x={x}")
        # decode & check
        _, _, _, lid_vars = vector_to_design(x)
        #print("DEBUG: LID vars:", lid_vars)
        
    
        
        
        # ---- LID test ----
        if not check_lid_constraints(lid_vars, catchment_areas):
            # print("LID constraints failed")
            lid_fail += 1
            continue

        # ---- connectivity (structure) test ----
        # ---- connectivity test on the enabled subgraph ----
        enabled_edges = [
            (u,v) for (u,v),idx in edge_to_pipe_idx.items() 
                    if pipe_enabled[idx]
        ]
        G_enabled = G_full.edge_subgraph(enabled_edges).to_undirected()

        # 1) 整体要连通
        if not nx.is_connected(G_enabled):
            #print("Warning: enabled pipes do not form a connected graph.")
            struct_fail += 1
            continue

        # 2) 验证每个出水口仍在子图中，且两两互联
        disconnected = False
        for o in sel:
            # 2.1) 节点本身要在启用子图中
            if not G_enabled.has_node(o):
                struct_fail += 1
                disconnected = True
                print("disconnected node")
                break
            # 2.2) 必须能与另一个出水口相互连通
            for t in sel:
                if o == t:
                    continue
                # 先确认 t 也是 G_enabled 的节点，再测路径
                if not G_enabled.has_node(t) or not nx.has_path(G_enabled, o, t):
                    struct_fail += 1
                    disconnected = True
                    print("break 2")
                    break
            if disconnected:
                print("break 3")
                break

        if disconnected:
            continue


        
        print(f"✓ Found valid design after {attempt+1} tries "
              f"(no_outlet={no_system_outfall}, lid_fail={lid_fail}, struct_fail={struct_fail})")
        # print the last x
        #print("DEBUG: Final x:", x)
        
        return x
    
# SWMM Simulation Execution Function
def generate_inp_file(template_inp_path, output_inp_path,
                      pipe_enabled, pipe_diameters, conduits_info,
                      outlet_selection, lid_vars,
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
    for row in conduits_mod:
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
    for row in conduits_mod:
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
    
def run_swmm_and_evaluate(inp_path, outfall_nodes, cost_per_diameter):
    

    # lcc
    lcc = compute_LCC(inp_path, cost_per_diameter, lid_vars,
                  cost_lid_bc=200, cost_lid_pp=150,
                  om_lid_bc=0.08, om_lid_pp=0.04,
                  discount_rate=0.02, project_lifetime=30)
    
    return lcc    

N_initial = 10
template_inp_path = "./AhvazNull.inp"
swmm_exe_path = "./swmm5.exe"  
max_iter = 40           # additional BO iterations (11th… up to 50th)
tol_ei   = 1e-3         # stop if max EI < tol

# loop to check if X_init is valid   
#X_init = np.array([sample_valid_design(conduits_info, catchment_areas, weight_list, weight_to_state) for _ in range(N_initial)], dtype=np.float32)
samples = []
for i in range(N_initial):
    X_init = sample_valid_design(conduits_info, catchment_areas, weight_list, weight_to_state) 
    if X_init is None or not isinstance(X_init, np.ndarray) or X_init.ndim != 1:
        print(f"Sample {i} returned bad x (None or not 1D), skipping:", X_init)
        continue
    if np.any(np.isnan(X_init)) or np.any(np.isinf(X_init)):
        print(f"Sample {i} contains NaN/Inf, skipping.")
        continue    
    samples.append(X_init)
if len(samples) == 0:
    raise RuntimeError("Failed to generate a valid initial design.")    
X_init = np.stack(samples, axis=0).astype(np.float32)
print(f"Collected {len(samples)} valid samples. X_init.shape = {X_init.shape}")

y_init = []
for i, x in enumerate(X_init):
    # if x in none or nan, print X_init
    #if np.any(np.isnan(x)) or np.any(np.isinf(x)):
    #    print(f"Warning: Sample {i} contains NaN or Inf values, skipping.")
        # print X_init
    #    print ("=============",X_init, "=============") 
    # Decode the design.
    pipe_enabled, pipe_diameters, outlet_selection, lid_vars = vector_to_design(x)
    #print("DEBUG lid_vars:", lid_vars, "shape:", lid_vars.shape)
    if not check_lid_constraints(lid_vars, catchment_areas):
        print("  → LID 约束不通过:", lid_vars[0::2] + lid_vars[1::2])
        continue

    conduits_mod = conduits_info
    coords_mod   = coords_info

    # 用 subcatch_outlets 作为新的 subcatchout_info
    #new_subcatchout_info = subcatch_outlets.copy()
    
    
    print(f">>> Running sample {i} …", flush=True)
    # Generate .inp file for the sample.
    out_inp = f"opt_trial_{i}.inp"
    rpt_file = f"opt_trial_{i}.rpt"
    out_file = f"opt_trial_{i}.out"
    
    try:
        
        generate_inp_file(template_inp_path, out_inp, pipe_enabled, pipe_diameters,
                            conduits_mod,        # 切割后
                            outlet_selection, lid_vars,loops, subcatch_info, coords_mod)
        
        # Run SWMM simulation.
        print("    -> About to call SWMM subprocess", flush=True)
        trunsim = time.time()
        lcc = run_swmm_and_evaluate(out_inp, outlet_selection, cost_per_diameter)
        print(f"    <- SWMM subprocess finished in {time.time()-trunsim:.1f}seconds")
        # Define objective as a combination. For example: y = qmax + alpha * lcc, alpha=1 for simplicity.
        qmax_list = []
        for i in range(N_initial):
            qmax = get_system_qmax_from_rpt(rpt_file)
            qmax_list.append(qmax)
        print(f"    <- SWMM returned qmax={qmax}, lcc={lcc}", flush=True)
        y = qmax + lcc
    except Exception as e:
        print("Simulation error on sample", i, ":", e)
        y = 1e6
    y_init.append(y)