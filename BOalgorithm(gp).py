import numpy as np
import pandas as pd
import random
import torch
import gpytorch
from gpytorch.priors import GammaPrior
from gpytorch.constraints import GreaterThan
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
import json
from concurrent.futures import ProcessPoolExecutor


# Optional: Random Forest surrogate
try:
    from sklearn.ensemble import RandomForestRegressor
except Exception:
    RandomForestRegressor = None

# Surrogate choice: 'gp' (default) or 'rf'
SURROGATE = os.environ.get('SURROGATE', 'gp').lower()

# Number of initial samples and BO iterations
N_initial = 200
max_iter = 200         

# Normalization bounds
LCC_MIN, LCC_MAX = 3e4, 3e6
Q_MAX_MIN, Q_MAX_MAX = 0, 3e4
EPS = 1e-6

# Objective weights (minimize w1*LCC_norm + w2*Qmax_norm)
w1 = 0.5  # weight for LCC
w2 = 0.5  # weight for max outflow

RAIN_TS_LIST = ["ssp585_Ens10y2h"]
'''
RAIN_TS_LIST = [
    "ACCESS-CM2_10y2h",
    "ACCESS-ESM1-5_10y2h",
    "BCC-CSM2-MR_10y2h",
    "CanESM5_10y2h",    # model 1
    "CESM2_10y2h",
    "CESM2-WACCM_10y2h",
    "CMCC-CM2-SR5_10y2h",
    "CMCC-ESM2_10y2h",
    "CNRM-CM6-1_10y2h",
    "CNRM-ESM2-1_10y2h",
    "EC-Earth3_10y2h",
    "EC-Earth3-Veg-LR_10y2h",
    "FGOALS-g3_10y2h",
    "GFDL-CM4_10y2h",
    "GFDL-ESM4_10y2h",
    "GISS-E2-1-G_10y2h",
    "HadGEM3-GC31-LL_10y2h",
    "HadGEM3-GC31-MM_10y2h",
    "IITM-ESM_10y2h",
    "INM-CM4-8_10y2h",
    "INM-CM5-0_10y2h",
    "IPSL-CM6A-LR_10y2h",
    "KACE-1-0-G_10y2h",
    "KIOST-ESM_10y2h",
    "MIROC6_10y2h",
    "MIROC-ES2L_10y2h",
    "MPI-ESM1-2-HR_10y2h",
    "MRI-ESM2-0_10y2h",
    "NESM3_10y2h",
    "NorESM2-MM_10y2h"  # model 2
]
'''

# Load conduits
conduits_info = np.loadtxt("[CONDUITS].txt")

# Load commercial diameters

commercial_df = pd.read_csv("Commercial.txt", sep=r"\s+", header=None)
commercial_df.columns = ["Diameter_m", "UnitCost"]
commercial_diameters = commercial_df["Diameter_m"].tolist()
# cost dictionary from file (if needed)
cost_per_diameter = dict(zip(commercial_df["Diameter_m"], commercial_df["UnitCost"]))

        

# Discrete optional weights (pipe diameters) list, strictly increasing
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
# each row: [subcatchment_ID, OutletNode]
subcatchout_info = np.loadtxt("Subcatchout.txt", dtype=int, delimiter=",")
# pipe_ids designated manually
pipe_ids = [2,17,12,10,8,15,16,23,24,30,32,37,44,55,41,48,53,50]
n_sc = subcatch_info.shape[0]
# create a dict: {ID → OutletNode}
#outlet_map = { row[0]: row[1] for row in subcatchout_info }
coords_info = np.loadtxt("[COORDINATES].txt", dtype=float, delimiter=",").reshape(-1, 3)

junctions_info = np.loadtxt("[JUNCTIONS].txt", dtype=float)
elev_map = { int(row[0]): float(row[1]) for row in junctions_info }

all_nodes = np.unique(conduits_info[:,1:3].astype(int)) # number of nodes: 48
# node_list: keep order stable
node_list = sorted(all_nodes)
node_to_idx = { node:i for i,node in enumerate(node_list) }
M = len(node_list)      # number of outfall: 48
LID_LEVELS = (0.0, 0.025, 0.05, 0.10)


def static_split(conduits_info, coords_info, outlet_nodes, pipe_ids):
    """
    Statically insert outlet_nodes[i] into pipe pipe_ids[i] within conduits_info,
    producing a master conduits_info and corresponding coords_info.
    Requires outlet_nodes and pipe_ids to have the same length.
    Returns master_coords_info and master_conduits_info.
    """
    coords_dict = {int(r[0]): (r[1], r[2]) for r in coords_info}
    master_nodes = coords_info.copy().tolist()
    master_conduits = [row.copy() for row in conduits_info]

    next_node_id = int(coords_info[:,0].max()) + 1
    next_pipe_id = int(conduits_info[:,0].max()) + 1

    for outlet_node, orig_pipe in zip(outlet_nodes, pipe_ids):
    # 1) find the conduit and its index in master_conduits
        for idx, row in enumerate(master_conduits):
            if int(row[0]) == orig_pipe:
                break
        else:
            raise RuntimeError(f"unable to find conduit ID={orig_pipe}")
        u, v = int(row[1]), int(row[2])
        xu, yu = coords_dict[u]; xv, yv = coords_dict[v]
        xo, yo = coords_dict[outlet_node]

    # 2) interpolate new node
    # if outlet_node in coords_dict, just use it without creating a new one
        if outlet_node in coords_dict:
            new_id = outlet_node
        else:
            new_id = next_node_id
            next_node_id += 1
            master_nodes.append([new_id, xo, yo])

    # 3) split the conduit u->v into u->new_id and new_id->v
    #    use original pipe id for the first segment, and new pipe id for the second segment
        row1 = row.copy()
        row1[2] = new_id  # u->new_id
        row2 = row.copy()
        row2[0] = next_pipe_id; next_pipe_id += 1
        row2[1] = new_id  # new_id->v

    # 4) use pop(idx) to remove the original row, then insert the two segments at the same position
        master_conduits.pop(idx)
        master_conduits.insert(idx, row1)
        master_conduits.insert(idx+1, row2)

    return np.array(master_nodes), np.array(master_conduits)

# 1) static split to obtain master network
master_coords, master_conduits = static_split(
    conduits_info,
    coords_info,
    outlet_nodes=subcatchout_info[:,1].tolist(),
    pipe_ids=pipe_ids
)

# 2) overwrite the original
coords_info = master_coords
conduits_info = master_conduits

# update global node dictionary
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


# 1) add every MST pipe-index as a node (index)
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

# subcatch_outlets
subcatch_outlets = set(subcatchout_info[:,1].tolist())

# find all nodes in G_full with out-degree 0
candidates = [n for n in G_full.nodes() if G_full.out_degree(n)==0]

# exclude subcatchment outlet nodes, keeping only true network endpoints
system_outfalls = [n for n in candidates if n not in subcatch_outlets]
if not system_outfalls:
    raise RuntimeError("No valid system outfall found!")
# randomly select two as outfalls
system_outfall_node = random.sample(system_outfalls, 2)

OUTLET_FIXED_MASK = np.zeros(M, dtype=np.uint8)
for n in system_outfall_node:
    OUTLET_FIXED_MASK[node_to_idx[n]] = 1
    
D = N_pipes + M + 36 # dimension of variables (dynamic)
#d_cont = 36
#d_disc = 65 + 18 + M
d_cont = N_pipes 
d_disc = M + 36

def eval_one_scenario_parallel(args):
    """
    run SWMM for a specific RAIN_TS scenario in parallel and return Qmax
    """
    if len(args) == 14:
        (ts_name,
         base_inp, ext_inp,
         base_rpt, ext_rpt,
         pipe_enabled,
         pipe_diameters,
         lid_vars,
         conduits_info,
         loops,
         subcatch_info,
         coords_info,
         template_inp_path,
         swmm_exe_path) = args
        run_uid = None
    elif len(args) == 15:
        (ts_name,
         base_inp, ext_inp,
         base_rpt, ext_rpt,
         pipe_enabled,
         pipe_diameters,
         lid_vars,
         conduits_info,
         loops,
         subcatch_info,
         coords_info,
         template_inp_path,
         swmm_exe_path,
         run_uid) = args
    else:
        raise ValueError(f"Unexpected args length={len(args)} in eval_one_scenario_parallel")

    # each scenario has its own independent file
    pid = os.getpid()
    uid = run_uid or f"pid{pid}_{int(time.time()*1000)}"
    inp_sc = f"{base_inp}_{ts_name}_{uid}{ext_inp}"
    rpt_sc = f"{base_rpt}_{ts_name}_{uid}{ext_rpt}"


    #inp_sc = f"{base_inp}_{ts_name}{ext_inp}"
    #rpt_sc = f"{base_rpt}_{ts_name}{ext_rpt}"

    # generate .inp
    generate_inp_file(
        template_inp_path,
        inp_sc,
        pipe_enabled,
        pipe_diameters,
        conduits_info,
        lid_vars,
        loops,
        subcatch_info,
        coords_info,
        rain_gage=ts_name,
    )

    # run SWMM
    subprocess.run(
        [swmm_exe_path, inp_sc, rpt_sc],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )

    # parse the Qmax for this scenario
    q_here = get_system_qmax_from_rpt(rpt_sc)
    return float(q_here)

def repair_monotone_global(pipe_enabled: np.ndarray,
                           pipe_diameters: List[Optional[float]],
                           weight_list: List[float],
                           edge_to_pipe_idx: Dict[Tuple[int, int], int],
                           n_iter: int = 10) -> None:
    """Global 'upstream ≤ downstream' repair: first map diameter to discrete code, then raise downstream in multiple rounds. In-place modification of pipe_diameters."""
    L = len(weight_list)
    code_of = {w: i for i, w in enumerate(weight_list)}
    codes = np.full(len(pipe_enabled), -1, dtype=int)
    for i, (en, d) in enumerate(zip(pipe_enabled, pipe_diameters)):
        if en and (d in code_of):
            codes[i] = code_of[d]
    for _ in range(n_iter):
        changed = False
        # notice: here assume that edge_to_pipe_idx contains all directed edges (u,v)
        for (u, v), idx_d in edge_to_pipe_idx.items():
            if not pipe_enabled[idx_d] or codes[idx_d] < 0:
                continue
            # all "direct upstream" edges (w,u)
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
    """Quantize to 4 levels and enforce BC+PP ≤ 0.10."""
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
            # keep the greater one, reduce the smaller one to the maximum level not exceeding the remaining quota
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
    fixed_outlet_mask: Optional[np.ndarray] = None,  # If you want to fix the system outlet like GA, you can pass a 0/1 mask of shape (M_nodes,)
    #mst_indices: Optional[Sequence[int]] = None,
)-> tuple[np.ndarray, list[Optional[float]], np.ndarray, np.ndarray]:
    """
    Decoder: BO continuous vectors x → (pipe_enabled, pipe_diameters, outlet_bits, lid_vars)
    x = [states_norm( N_pipes ), outlet_bits( M_nodes ), lid_vars( 2*S )], where states_norm ∈ [0,1] represents "diameter level / L", not 0/1.
    """
    L = len(weight_list)
    assert x.shape[0] >= N_pipes + M_nodes, "x dimension is insufficient"

    # 1) diameter level decoding (note: no longer using *10 hardcoding)
    states_norm = np.clip(x[:N_pipes], 0.0, 1.0)
    states = np.clip(np.round(states_norm * L), 0, L).astype(int)  # 0..L
    pipe_enabled = (states > 0).astype(int)

    pipe_diameters: List[Optional[float]] = [None] * N_pipes
    for i, s in enumerate(states):
        if s >= 1:
            pipe_diameters[i] = weight_list[s - 1]

    # 2) outlet bits
    outlet_slice = x[N_pipes : N_pipes + M_nodes]
    outlet_bits = np.round(np.clip(outlet_slice, 0.0, 1.0)).astype(int)
    if fixed_outlet_mask is not None:
    # Align with GA: if the system outlet must be 1, then cover or take the maximum
        outlet_bits = np.maximum(outlet_bits, fixed_outlet_mask.astype(int))
        # If strictly only allowing fixed to be 1, you can change to: outlet_bits = fixed_outlet_mask.astype(int).copy()

    # 3) LID variables (BC, PP, BC, PP, ...)
    lid_vars = x[N_pipes + M_nodes : ].astype(float).copy()
    lid_vars = quantize_lid_vars_4levels(lid_vars)

    # 4) (optional) global monotonicity repair: upstream ≤ downstream
    if enforce_monotone and edge_to_pipe_idx is not None:
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
                             require_single_component=True,     # the entire network is a weakly connected component
                             require_reachability=False          # each node can reach at least one outlet
                             ):

    if not outlet_nodes:                  # at least one outlet must exist
        return False

    # Enabled subgraph (directed)
    enabled_edges = [(u,v) for (u,v),idx in edge_to_pipe_idx.items() if pipe_enabled[idx]]
    if not enabled_edges:
        return False
    H = G_full.edge_subgraph(enabled_edges).copy()
    U = H.to_undirected()

    # Outlets must be in the subgraph, and must be sinks (out-degree=0) in H, and not belong to subcatchment interfaces
    for o in outlet_nodes:
        if o not in H: return False
        if H.out_degree(o) != 0: return False
        if subcatch_outlets and o in subcatch_outlets: return False

    # Weak connectivity: the entire network is a single component, or relaxed to "each component has at least one outlet"
    if require_single_component:
        if not nx.is_connected(U):
            return False
    else:
        for comp_nodes in nx.connected_components(U):
            if not any(o in comp_nodes for o in outlet_nodes):
                return False

    # Upstream ≤ Downstream (only for enabled edges)
    for (u,v), idx_d in edge_to_pipe_idx.items():
        if not pipe_enabled[idx_d]:
            continue
        d_down = pipe_diameters[idx_d]
        for (w,uu), idx_u in edge_to_pipe_idx.items():
            if uu == u and pipe_enabled[idx_u]:
                d_up = pipe_diameters[idx_u]
                if d_up > 0 and d_down > 0 and d_up > d_down + 1e-9:
                    return False

    # Reachability check: perform multi-source search in the "reverse graph" from all outlets
    if require_reachability:
        Hr = H.reverse(copy=False)
        seen = set(outlet_nodes)
        for s in outlet_nodes:
            seen |= nx.descendants(Hr, s) | {s}   # points that can reach s (can flow to s in the original graph)
        # only require reachability for nodes that appear in enabled edges
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

def debug_print_raingages(inp_path: str):
    """Print the [RAINGAGES] section of an .inp file (for debugging)."""
    print(f"\n--- [RAINGAGES] in {inp_path} ---")
    in_rg = False
    with open(inp_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not in_rg:
                # Match [RAINGAGE] or [RAINGAGES], case-insensitive
                if line.upper().startswith("[RAINGAGE"):
                    in_rg = True
                    print(raw.rstrip())
                continue
            # if we hit another section, stop
            if line.startswith("["):
                break
            print(raw.rstrip())
    print("--- end [RAINGAGES] ---\n")

def init_population_feasible_bo(
    N_pipes: int,
    M_nodes: int,
    mst_indices: List[int],                             # MST to pipe index
    edge_to_pipe_idx: Dict[Tuple[int, int], int],       # (u,v)->idx
    node_to_idx: Dict[int, int],                        # node ID to 0..M-1 (for outlet bits)
    system_outfall_nodes: List[int],                    # system outfall node list (fixed as outlets)
    weight_list: List[float],                           # optional diameter grades (sorted from small to large)
    catchment_areas: List[float],                       # subcatchment areas (length = number of subcatchments S)
    add_prob: float = 0.15,                             # non-MST random edge addition probability (align with GA)
    rng: Optional[np.random.Generator] = None,
    max_try: int = 200,
    ensure_incoming_to_outfalls: bool = True   # if construction_check needed, ensure at least one incoming edge for each outlet
):
        """
        Generate a feasible individual in one step, equivalent to GA's init_population_feasible:
            - Enable edges: MST must be on + non-skeleton edges opened with probability add_prob
            - Diameters: assign random levels to enabled edges, then perform global monotone repair (upstream ≤ downstream)
            - Outlets: fix system endpoints to 1, others to 0
            - LID: generate (BC, PP) per subcatchment and ensure BC+PP ≤ 0.10
        Returns:
            pipe_enabled: (N_pipes,) int {0,1}
            pipe_diameters: List[float or None]
            outlet_bits: (M_nodes,) int {0,1}
            lid_vars: (2*S,) float  [BC1,PP1, BC2,PP2, ...]
            states_norm: (N_pipes,) float  # level normalized to [0,1] for BO vectorization
            x: (N_pipes + M_nodes + 2*S,)  # design vector (consistent with vector_to_design)
        """
    if rng is None:
        rng = np.random.default_rng()
    L = len(weight_list)

    for _ in range(max_try):
        
    # 1) Enable edges: MST must be on + non-skeleton edges opened with add_prob
        pipe_enabled = np.zeros(N_pipes, dtype=int)
        if len(mst_indices) > 0:
            pipe_enabled[np.array(mst_indices, dtype=int)] = 1
        flip = rng.random(N_pipes) < add_prob
        pipe_enabled = np.where((pipe_enabled == 1) | flip, 1, 0).astype(int)
        '''
        # (if only MST): 确保每个系统末端至少一条入边启用（若骨架未覆盖到末端，启用最近一条）
        if ensure_incoming_to_outfalls:
            for v in system_outfall_nodes:
                incoming = [idx for (u, w), idx in edge_to_pipe_idx.items() if w == v]
                if incoming and not any(pipe_enabled[idx] == 1 for idx in incoming):
                    pipe_enabled[incoming[0]] = 1
        '''
    # 2) Initial diameter level: for enabled edges sample in [1..L], disabled edges use level 0 (placeholder)
        '''
        states = np.where(pipe_enabled == 1, L, 0).astype(int)  # 1..L；最大档= L
        pipe_diameters = [weight_list[s-1] if s > 0 else None for s in states]
        states_norm = states / float(L)  # 归一化给 BO 用；注意后面会再修复，但归一化的向量也要同步更新
        '''
        states = np.zeros(N_pipes, dtype=int)
        k = np.where(pipe_enabled == 1)[0]
        if k.size > 0:
            states[k] = rng.integers(low=1, high=L + 1, size=k.size, endpoint=False) + 0  # [1..L]
    # Normalize to [0,1] for BO; later we repair but keep normalized vector in sync
        states_norm = states / float(L)
        
    # 3) Outlet bits: fix system endpoints to 1, others to 0
        '''
        outlet_bits = np.zeros(M_nodes, dtype=int)
        for n in system_outfall_nodes:
            j = node_to_idx.get(n, None)
            if j is not None: outlet_bits[j] = 1

        '''
        if OUTLET_FIXED_MASK is not None:
            outlet_bits = OUTLET_FIXED_MASK.astype(int).copy()
            assert outlet_bits.shape == (M_nodes,)
        else:
            outlet_bits = np.zeros(M_nodes, dtype=int)
            for n in system_outfall_nodes:
                j = node_to_idx.get(n, None)
                if j is not None: outlet_bits[j] = 1
        
    # 4) LID: discrete
        lid_pairs = []
        '''
        sample_pair = _sample_lid_pair_continuous 
        for _ in range(len(catchment_areas)):
            lid_pairs.append(sample_pair())
    # Flatten to [BC1, PP1, BC2, PP2, ...]
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
        
    # 5) Map states → diameters and perform global monotone repair
        pipe_diameters: List[Optional[float]] = [None] * N_pipes
        for i, s in enumerate(states):
            if s >= 1:
                pipe_diameters[i] = weight_list[s - 1]
        repair_monotone_global(pipe_enabled, pipe_diameters, weight_list, edge_to_pipe_idx, n_iter=10)

    # 6) Build BO design vector x (note: diameter levels are ratios, not 0/1)
        x = np.concatenate([states_norm.astype(np.float32),
                            outlet_bits.astype(np.float32),
                            lid_vars.astype(np.float32)], axis=0)
        
    # Structural constraints: support bool or (bool, msg)
        pipe_enabled, pipe_diams, outlet_bits, lid_vars = vector_to_design(x,
                                                                            N_pipes=N_pipes,
                                                                            M_nodes=M,
                                                                            weight_list=weight_list,
                                                                            enforce_monotone=True,
                                                                            edge_to_pipe_idx=edge_to_pipe_idx,
                                                                            fixed_outlet_mask=OUTLET_FIXED_MASK,  # key
                                                                            #mst_indices=mst_indices
                                                                            )
        idx_to_node = {v: k for k, v in node_to_idx.items()}
        outlet_idx = np.flatnonzero(np.asarray(outlet_bits)).tolist()
    outlet_nodes_list = [idx_to_node[j] for j in outlet_idx]  # Python list

        _struct_ret = check_structure_validity(pipe_enabled, pipe_diams, outlet_nodes_list,
                                               edge_to_pipe_idx, G_full,
                                               subcatch_outlets=None,
                                               require_single_component=True,     # entire network is one weakly connected component
                                               require_reachability=False          # each node can reach at least one outlet
                                               )
        if isinstance(_struct_ret, tuple):
            ok_struct = bool(_struct_ret[0])
            msg_s = _struct_ret[1] if len(_struct_ret) > 1 else ""
        else:
            ok_struct = bool(_struct_ret)
            msg_s = ""

    # LID constraints: support different signatures & return forms
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
                      rain_gage=None
                      ):
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
    # Ensure: each node in outfall_selection is pointed to by at most one conduit; and never appears as from_node
    used_outfalls = set() 
    #emitted = []              # <— keep track of which pipe_ids were actually written
    for row in conduits_info:
        pipe_id = int(row[0])
        frm = int(row[1]); to = int(row[2])
    # 1) Keep only allowed conduits
        idx = edge_to_pipe_idx.get((frm, to))
        if idx is None or not pipe_enabled[idx]:
            continue
    # 2) If this conduit’s to_node is a system outfall and already used once, skip it
        if to in system_outfall_node:
            if to in used_outfalls:
                # Disable any second or more conduits to the same outfall node
                continue
            # First encounter, record
            used_outfalls.add(to)

    # 3) Also ensure this outfall will not be used as from_node
        if frm in system_outfall_node:
            # If it is an outfall, it must never be the start of a conduit
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
    # Find global idx, check if it’s enabled
        idx = edge_to_pipe_idx.get((frm, to))
        if idx is None or not pipe_enabled[idx]:
            continue

    # Likewise: an outfall can only be used once
        if to in system_outfall_node:
            if to in used_outfalls:
                continue
            used_outfalls.add(to)

    # Outfall cannot have outgoing degree
        if frm in system_outfall_node:
            continue
        d = pipe_diameters[idx]
        
        '''
    # If executable, comment out
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
    sid     = i+67         # subcatchment ID, e.g., 1 → S1
    outlet  = int(subcatchout_info[i, 1])         # template’s expected outlet node ID
        row      = subcatch_info[i]
    rg      = int(row[1])         # rain gage ID
    area    = row[3]              # area
    imperv  = row[4]              # imperviousness (%)
    width   = row[5]              # width
    slope   = row[6]              # slope (%)
    curb    = row[7]              # curb length
    #snow    = row[8]              # snow pack (mm)
    # SnowPack filled with 0
        subcat.append(
            f"{sid:<14d}"   # 左对齐，占14列
            f"{rg:<16d}"    # 左对齐，占16列
            f"{outlet:<16d}"# 左对齐，占16列
            f"{area:>8.2f}" # 右对齐，占8列，保留2位小数
            f"{imperv:>8.2f}"
            f"{width:>8.2f}"
            f"{slope:>8.2f}"
            f"{curb:>8.2f}"
            f"\n"            # empty SnowPack column
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
    # sid is the new ID (67…84) for column 1
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
    unit_area_bc  = 200     # theoretical max area per BC unit
    unit_area_pp  = 150 
    fromImp_BC = 57.14286
    fromImp_PP = 42.85714
    bc_fracs = lid_vars[0::2]   # length 18, fraction of each catchment’s area used by BC
    pp_fracs = lid_vars[1::2]   # length 18, fraction for PP
    for i, (bc_frac, pp_frac, A) in enumerate(zip(bc_fracs, pp_fracs, catchment_areas)):
        
        sid = start_id + i
    # --- BC unit ---
    # Assume each subcatchment has at most one BC unit (Number=1),
    # If you want dynamic num_bc, read lv["num_bc"] from lid variables
        bc_area = bc_frac * A
        if bc_area > 0:
            num_bc       = int(np.ceil(bc_area / unit_area_bc)) if bc_area > 0 else 0
            #area_per_bc  = bc_area / num_bc
            from_imp_bc  = fromImp_BC if num_bc > 0 else 0.0
            # Column 7 FromImp = unit area / subcatchment area * 100 (%)
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
    # --- PP unit ---
        pp_area = pp_frac * A
        if pp_area > 0:
            num_pp       = int(np.ceil(pp_area / unit_area_pp)) if pp_area > 0 else 0
            #area_per_pp  = pp_area / num_pp
            from_imp_pp  = fromImp_PP if num_pp > 0 else 0
            lidu.append(
                f"{sid:<14d}"
                f"{2:<16d}"        # LID Process: assign code 2 for PP
                f"{num_pp:<9d}"
                f"{unit_area_pp:>10.2f}"
                f"{width:>11d}"
                f"{0:>11d}"
                f"{from_imp_pp:>11.5f}"
                f"{1:>11d}"
                "\n"
            )
    #print("DEBUG lid_vars:", lid_vars)      
    # Write back LID_USAGE block
    lines = replace_block(lines, "LID_USAGE", lidu)
  
    # 9) Junctions
    
    junc_block = [
        ";;Name           Elevation  MaxDepth   InitDepth  SurDepth   Aponded   \n",
        ";;-------------- ---------- ---------- ---------- ---------- ----------\n"
    ]
    for row in junctions_info:
        node = int(row[0])
        if node in  system_outfall_node:
            # If an outfall node, skip
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
            f"{int(node):<16d}"   # Node ID, left-aligned in 16 chars
            f"{x:>12.1f}"         # X-Coord, right-aligned, 1 decimal
            f"{y:>13.1f}"         # Y-Coord, right-aligned, 1 decimal
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
    
    # Require explicit rain_gage; caller should loop RAIN_TS_LIST and pass each name
    if rain_gage is None:
        raise ValueError("generate_inp_file requires a rain_gage name; pass each RAIN_TS_LIST item explicitly.")

    # Update [RAINGAGES] to reference TIMESERIES <rain_gage>
    rg_start = None
    for i, raw in enumerate(lines):
        hdr = raw.strip().upper()
        if hdr.startswith("[RAINGAGE"):
            rg_start = i + 1
            break
    if rg_start is not None:
        j = rg_start
        while j < len(lines) and not lines[j].strip().startswith("["):
            raw = lines[j]
            stripped = raw.strip()
            if stripped and not stripped.startswith(";"):
                if "TIMESERIES" in stripped.upper():
                    lines[j] = re.sub(r'(?i)(TIMESERIES\s+)\S+', r'\1' + str(rain_gage), raw)
            j += 1
    '''
    # Ensure [TIMESERIES] has an entry named rain_gage; rename the first data line
    ts_start = None
    for i, raw in enumerate(lines):
        if raw.strip().upper().startswith("[TIMESERIES]"):
            ts_start = i + 1
            break
    if ts_start is None:
        raise RuntimeError("[TIMESERIES] section not found; cannot reference undefined series from [RAINGAGES].")
    j = ts_start
    renamed = False
    while j < len(lines) and not lines[j].strip().startswith("["):
        raw = lines[j]
        stripped = raw.strip()
        if stripped and not stripped.startswith(";;") and not stripped.startswith(";"):
            # Preserve leading whitespace and the remainder of the line
            m = re.match(r"^(\s*)", raw)
            leading = m.group(1) if m else ""
            tokens = re.split(r"\s+", stripped)
            if tokens:
                name_len = len(tokens[0])
                rest = raw[len(leading) + name_len:]
                lines[j] = f"{leading}{rain_gage}{rest}"
                renamed = True
                break
        j += 1
    if not renamed:
        raise RuntimeError("[TIMESERIES] section has no data lines to rename; add entries or provide format.")
    '''

    # write out
    with open(output_inp_path, 'w') as f:
        f.writelines(lines)

    return output_inp_path

def get_system_qmax_from_rpt(rpt_path):
    """Extract System max flow from an existing rpt file without running SWMM."""
    with open(rpt_path) as f:
        for raw in f:
            line = raw.strip()
            # Case-insensitive: if line starts with System, parse the 4th column
            if re.match(r'(?i)^system\b', line):
                parts = re.split(r'\s+', line)
                if len(parts) >= 4:
                    try:
                        return float(parts[3])
                    except ValueError:
                        pass
    raise RuntimeError(f"System max flow not found in {rpt_path}.")

def compute_LCC(inp_path,cost_per_diameter,lid_vars,cost_lid_bc=200,cost_lid_pp=150,om_lid_bc=0.08,om_lid_pp=0.04,discount_rate=0.02,project_lifetime=30):
     """
     1) Read [XSECTIONS]: 3rd column is diameter → dict name->diam
     2) Read [CONDUITS]: 4th column is Length, 1st column is Link Name
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
                    # parts[2] is the 3rd column diameter
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
                    length = float(parts[3])   # 4th column, Length
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
        Feasible region projection for continuous encoding:
            - Raise states_norm for must-on edges (MST, connectivity) to >= min_level/L
            - For robustness, optionally set these edges to 'max' (=1.0), analogous to old "1→10 levels"
            - Ensure 2 outlets from system_outfall_nodes
            - LID: apply upper-cap scaling or 4-level quantization + repair.
        Returns x_fixed.
    """
    x = np.array(x_raw, dtype=float).copy()
    assert x.ndim == 1

    # 1) 切段
    states_norm = x[:N_pipes]                    # 连续档位比例
    out_bits    = x[N_pipes:N_pipes+M_nodes]
    lid_vars    = x[N_pipes+M_nodes:]

    # 2) Outlets: ensure two nodes from outfall set
    chosen = [node_list[i] for i,b in enumerate(out_bits>0.5) if b]
    # 补足/裁剪
    chosen = list(dict.fromkeys([n for n in chosen if n in system_outfall_nodes]))[:2]
    import random
    while len(chosen) < 2:
        c = random.choice(system_outfall_nodes)
        if c not in chosen:
            chosen.append(c)
    out_bits_fixed = np.array([1 if n in chosen else 0 for n in node_list], dtype=float)

    # 3) Must-on edges: MST + paths ensuring each outlet connects to chosen[0]
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

    # 4) Raise states_norm of these must-on edges above 0
    if force_level == 'max':
        target = 1.0
    elif force_level == 'mid':
        target = (min_level + L) / (2*L)  # 粗略中位比例
    else:  # 'min'
        target = max(min_level / L, 1e-6)

    states_norm_fixed = states_norm.copy()
    for idx in must_on:
        states_norm_fixed[idx] = max(states_norm_fixed[idx], target)

    # 5) LID: apply cap scaling (or quantize to 4 levels then repair ≤0.10)
    # 5.1 Direct scaling version:
    lid = lid_vars.copy()
    for i in range(len(catchment_areas)):
        bc, pp = lid[2*i], lid[2*i+1]
        tot = bc + pp
        if tot > max_lid_frac:
            fac = max_lid_frac / tot
            lid[2*i] *= fac; lid[2*i+1] *= fac

    # If already converted to four levels, you can replace with:
    # lid = quantize_lid_vars_4levels(lid)  # includes ≤0.10 repair

    # 6) 写回
    x_fixed = x.copy()
    x_fixed[:N_pipes] = np.clip(states_norm_fixed, 0.0, 1.0)
    x_fixed[N_pipes:N_pipes+M_nodes] = out_bits_fixed
    x_fixed[N_pipes+M_nodes:] = lid
    return x_fixed


twhole = time.time()  # [TIMING] total time    
class MixedKernel(gpytorch.kernels.Kernel):
    """Original multiplicative kernel: k = k_cont * k_disc.
    Kept for backward compatibility / ablation.
    """
    def __init__(self, cont_kernel, disc_kernel):
        super().__init__()
        self.cont_kernel = cont_kernel
        self.disc_kernel = disc_kernel

    def forward(self, x1, x2, diag=False, **params):
        x1_cont = x1[..., :d_cont]; x2_cont = x2[..., :d_cont]
        x1_disc = x1[..., d_cont:]; x2_disc = x2[..., d_cont:]
        k_cont = self.cont_kernel(x1_cont, x2_cont, diag=diag, **params)
        k_disc = self.disc_kernel(x1_disc, x2_disc, diag=diag, **params)
        return k_cont * k_disc

class AdditiveMixedKernel(gpytorch.kernels.Kernel):
    """Additive + optional interaction kernel.

    k(x,x') = k_cont(x_c,x'_c) + k_disc(x_d,x'_d) (+ interaction * weight)
    The interaction term helps restore some coupling lost when switching
    from a pure product to pure sum. Lightweight and more flexible than
    strict product in high sparse binary dimensions.
    """
    def __init__(self, cont_kernel, disc_kernel, include_interaction=True, interaction_weight=1.0):
        super().__init__()
        self.cont_kernel = cont_kernel
        self.disc_kernel = disc_kernel
        self.include_interaction = include_interaction
        self.register_parameter(name="raw_interaction_weight", parameter=torch.nn.Parameter(torch.tensor(float(interaction_weight))))

    def forward(self, x1, x2, diag=False, **params):
        x1_cont = x1[..., :d_cont]; x2_cont = x2[..., :d_cont]
        x1_disc = x1[..., d_cont:]; x2_disc = x2[..., d_cont:]
        k_cont = self.cont_kernel(x1_cont, x2_cont, diag=diag, **params)
        k_disc = self.disc_kernel(x1_disc, x2_disc, diag=diag, **params)
        k_sum = k_cont + k_disc
        if self.include_interaction:
            # softplus to keep weight positive, prevents sign flips harming PSD
            interaction_w = torch.nn.functional.softplus(self.raw_interaction_weight)
            return k_sum + interaction_w * (k_cont * k_disc)
        return k_sum

# Flag to choose kernel style; set to 'add' to enable additive variant.
KERNEL_MODE = 'add'  # options: 'prod' (original), 'add'

# Sub-kernels:
cont_kernel = gpytorch.kernels.RBFKernel(ard_num_dims=d_cont)
disc_kernel = gpytorch.kernels.RBFKernel(ard_num_dims=d_disc)

# Initialize discrete lengthscale larger (avoid immediate near-diagonal kernel)
with torch.no_grad():
    disc_kernel.lengthscale.fill_(2.0)

# Register lengthscale priors (Gamma(k, theta=1/rate) style in gpytorch: Gamma(shape, rate))
cont_kernel.register_prior("ls_cont_prior", GammaPrior(2.0, 1.0), "lengthscale")
disc_kernel.register_prior("ls_disc_prior", GammaPrior(2.0, 1.0), "lengthscale")

if KERNEL_MODE == 'add':
    # Expose inner ScaleKernels to attach outputscale priors
    cont_scale_kernel = gpytorch.kernels.ScaleKernel(cont_kernel)
    disc_scale_kernel = gpytorch.kernels.ScaleKernel(disc_kernel)
    cont_scale_kernel.register_prior("os_cont_prior", GammaPrior(1.5, 1.0), "outputscale")
    disc_scale_kernel.register_prior("os_disc_prior", GammaPrior(1.5, 1.0), "outputscale")
    base_mixed = AdditiveMixedKernel(
        cont_scale_kernel,
        disc_scale_kernel,
        include_interaction=True,
        interaction_weight=0.5,
    )
else:
    base_mixed = MixedKernel(cont_kernel, disc_kernel)

combined_kernel = gpytorch.kernels.ScaleKernel(base_mixed)
combined_kernel.register_prior("os_global_prior", GammaPrior(1.5, 1.0), "outputscale")

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

    # Continuous: first N_pipes
    states_norm = X[:, :N_pipes]

    # Discrete: binarize outlet bits
    outlet_bits = (X[:, N_pipes:N_pipes+M_nodes] > 0.5).astype(np.float32)

    # Discrete: LID quantized to 4 levels (via unified function)
    lid = X[:, N_pipes+M_nodes : N_pipes+M_nodes + 2*S]
    lid_q = np.vstack([quantize_lid_vars_4levels(row) for row in lid])

    X_proj = np.concatenate([states_norm, outlet_bits, lid_q], axis=1)
    return X_proj[0] if single else X_proj


def reorder_batch(X):
    return np.array([reorder_for_kernel(x) for x in X])

def dedupe_xy(X: np.ndarray, y: np.ndarray):
    """Remove exact duplicate design rows in X; aggregate y by mean.
    Returns X_unique, y_agg with matching order.
    """
    if not isinstance(X, np.ndarray) or X.ndim != 2 or len(X) != len(y):
        return X, y
    # Use rounding to avoid float artifacts from projection
    Xr = np.round(X, 6)
    uniq, inv, counts = np.unique(Xr, axis=0, return_inverse=True, return_counts=True)
    if uniq.shape[0] == X.shape[0]:
        return X, y
    y_agg = np.zeros(uniq.shape[0], dtype=np.float32)
    for i in range(uniq.shape[0]):
        y_agg[i] = float(np.mean(y[inv == i]))
    return uniq.astype(np.float32), y_agg.astype(np.float32)

# --- Initial design utilities: Sobol sampling + maximin selection ---
def _sobol_pool(n_pool: int, dim: int, seed: int | None = None) -> np.ndarray:
    """Generate a Sobol pool in [0,1]^dim as float32 numpy array."""
    try:
        import torch
        from torch.quasirandom import SobolEngine
        eng = SobolEngine(dim, scramble=True, seed=seed)
        X = eng.draw(n_pool).cpu().numpy().astype(np.float32)
    except Exception:
        # Fallback to uniform if Sobol unavailable
        rng = np.random.default_rng(seed)
        X = rng.random((n_pool, dim), dtype=np.float32)
    return X

def _pairwise_min_dists(emb: np.ndarray, chosen_idx: list[int]) -> np.ndarray:
    """Compute min distance from each point to the chosen set (Euclidean)."""
    n = emb.shape[0]
    if not chosen_idx:
        return np.full(n, np.inf, dtype=np.float32)
    chosen = emb[chosen_idx]
    # Compute distances in a vectorized way
    dists = np.linalg.norm(emb[:, None, :] - chosen[None, :, :], axis=2)  # (n, k)
    return dists.min(axis=1)

def select_maximin(X_pool: np.ndarray, k: int) -> np.ndarray:
    """Greedy maximin selection of k points from X_pool.
    Distance is computed in the reordered/kernel space for better coverage of mixed variables.
    Returns indices of selected points (length k).
    """
    if X_pool.shape[0] <= k:
        return np.arange(X_pool.shape[0], dtype=int)
    emb = reorder_batch(X_pool).astype(np.float32)
    # Seed: farthest from the centroid
    centroid = emb.mean(axis=0, keepdims=True)
    d0 = np.linalg.norm(emb - centroid, axis=1)
    first = int(np.argmax(d0))
    chosen = [first]
    min_d = _pairwise_min_dists(emb, chosen)
    selected = set(chosen)
    while len(chosen) < k:
        # Exclude already selected by setting their score to -inf
        scores = min_d.copy()
        for idx in selected:
            scores[idx] = -np.inf
        nxt = int(np.argmax(scores))
        chosen.append(nxt)
        selected.add(nxt)
        # Update min distances incrementally
        d_new = np.linalg.norm(emb - emb[nxt], axis=1)
        min_d = np.minimum(min_d, d_new)
    return np.array(chosen, dtype=int)

# Define the Expected Improvement (EI) acquisition function.

def expected_improvement(X, f_best, model, xi=0.0):
    """
    X: (n, D) torch.float, already reordered
    f_best: float, current minimal objective (minimization)
    model: trained GP
    xi: exploration parameter (commonly 0~0.05; larger means more exploration)
    """
    model.eval()
    with torch.no_grad():
        post = model(X)
        mu   = post.mean            # (n,)
        var  = post.variance.clamp_min(1e-12)
        sigma = var.sqrt()

    # Minimization: imp = f_best - mu - xi
        imp = f_best - mu - xi
        Z   = imp / sigma

    # EI = imp * Phi(Z) + sigma * phi(Z); set EI=0 when sigma≈0
        normal = torch.distributions.Normal(0.0, 1.0)
        ei  = imp * normal.cdf(Z) + sigma * torch.exp(-0.5 * Z**2) / np.sqrt(2*np.pi)
        ei = torch.where(sigma <= 1e-12, torch.zeros_like(ei), ei)
    return ei

def train_rf_regressor(X_np: np.ndarray, y_np: np.ndarray):
    """Train a RandomForestRegressor on reordered features; returns a dict.
    Uses y z-score for stability; caller should compute acquisition on RF predictions.
    """
    if RandomForestRegressor is None:
        raise RuntimeError("sklearn not available: cannot use SURROGATE='rf'")
    X_feat = reorder_batch(X_np).astype(np.float32)
    y_mean = float(y_np.mean()); y_std = float(max(y_np.std(), 1e-6))
    y_z = (y_np - y_mean) / y_std
    rf = RandomForestRegressor(
        n_estimators=400,
        max_depth=None,
        min_samples_leaf=3,
        n_jobs=-1,
        random_state=42,
        oob_score=False,
    )
    rf.fit(X_feat, y_z)
    return {"rf": rf, "y_mean": y_mean, "y_std": y_std}

# If you want a piecewise schedule for xi in BO; comment out if not needed
def xi_schedule_piecewise(it, N_initial, max_iter,
                          xi_hi=0.1, xi_mid=0.05, xi_lo=0.01,
                          p1=0.33, p2=0.67):
    """it is the current iteration index (starting at N_initial); returns xi for this round."""
    k = max(0, it - N_initial)
    p = k / max(1, max_iter)  # 进度 0~1
    if p < p1:
        return xi_hi
    elif p < p2:
        return xi_mid
    else:
        return xi_lo
    


def Bayesian_Opt(inp_file, rpt_file, seed: int | None = None, plot: bool = True, w1_override=None, w2_override=None, tag: str = ""):
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)
        print(f"[BO] Seed = {seed}")
        # --- weights override (for sweep runs) ---
    w1_local = float(w1_override) if w1_override is not None else float(w1)
    w2_local = float(w2_override) if w2_override is not None else float(w2)

    # safety: normalize if user passes arbitrary numbers
    s = w1_local + w2_local
    if s <= 0:
        w1_local, w2_local = 0.5, 0.5
    else:
        w1_local /= s
        w2_local /= s

    if tag:
        print(f"[BO] tag={tag}  weights: w1(cost)={w1_local:.3f}, w2(flow)={w2_local:.3f}")

    samples = []
    q_init = []
    c_init = []
    tgenersample = time.time()  # [TIMING] sampling

    # 1) Generate a Sobol pool, project to feasible, dedupe
    n_pool = max(N_initial * 6, 240)
    tries = 0
    X_pool = None
    while tries < 3:
        tries += 1
        X_raw = _sobol_pool(n_pool, D, seed=12345 + tries)
        # Project each to feasible space
        X_proj = []
        for xr in X_raw:
            xp = project_to_feasible_continuous(
                xr,
                N_pipes=N_pipes, M_nodes=M, L=len(weight_list),
                node_list=node_list, conduits_info=conduits_info,
                mst_indices=sorted(mst_indices),
                system_outfall_nodes=system_outfall_node,
                catchment_areas=catchment_areas,
                max_lid_frac=0.10,
                force_level='mid'
            )
            if xp is not None and isinstance(xp, np.ndarray) and xp.ndim == 1 and not (
                np.any(np.isnan(xp)) or np.any(np.isinf(xp))
            ):
                X_proj.append(xp.astype(np.float32))
        if not X_proj:
            n_pool *= 2
            continue
        X_pool = np.vstack(X_proj)
    # Deduplicate by rounding
        Xr = np.round(X_pool, 6)
        _, uniq_idx = np.unique(Xr, axis=0, return_index=True)
        X_pool = X_pool[np.sort(uniq_idx)]
        if X_pool.shape[0] >= max(N_initial, 4):
            break
        n_pool *= 2

    if X_pool is None or X_pool.shape[0] == 0:
        raise RuntimeError("Failed to generate a feasible Sobol pool for initial design.")

    # 2) Maximin select N_initial designs
    sel_idx = select_maximin(X_pool, N_initial)
    X_init = X_pool[sel_idx]
    print(
        f"Generated {X_init.shape[0]} initial samples via Sobol+maximin "
        f"in {time.time() - tgenersample:.2f} seconds (pool={X_pool.shape[0]})."
    )
    idx = np.random.permutation(len(X_init))  # small shuffle for robustness
    tsampcheck = time.time()  # [TIMING] sample check

    # 3) Evaluate initial samples (compute LCC and SWMM Qmax)
    #    If RAIN_TS_LIST is non-empty → scenario-averaged Qmax over multiple TIMESERIES
    for i, X_i in enumerate(X_init):
        pipe_enabled, pipe_diams, outlet_bits, lid_vars = vector_to_design(
            X_i,
            N_pipes=N_pipes,
            M_nodes=M,
            weight_list=weight_list,
            enforce_monotone=True,
            edge_to_pipe_idx=edge_to_pipe_idx,
            fixed_outlet_mask=OUTLET_FIXED_MASK,
        )
        samples.append(X_i)

        if RAIN_TS_LIST:
            q_s_list = []
            lcc_val = None

            # derive base names so each scenario writes a different inp/rpt
            base_inp, ext_inp = os.path.splitext(inp_file)
            base_rpt, ext_rpt = os.path.splitext(rpt_file)

            '''
            for ts_name in RAIN_TS_LIST:
                # scenario-specific filenames
                inp_sc = f"{base_inp}_{ts_name}{ext_inp}"
                rpt_sc = f"{base_rpt}_{ts_name}{ext_rpt}"

                # Build .inp for this scenario: same design, different rainfall
                generate_inp_file(
                    template_inp_path,
                    inp_sc,
                    pipe_enabled,
                    pipe_diams,
                    conduits_info,
                    lid_vars,
                    loops,
                    subcatch_info,
                    coords_info,
                    rain_gage=ts_name,  # choose TIMESERIES name
                )
                debug_print_raingages(inp_sc)

                # LCC does not depend on rainfall → compute once
                if lcc_val is None:
                    lcc_val = compute_LCC(
                        inp_sc,
                        cost_per_diameter,
                        lid_vars,
                        cost_lid_bc=200,
                        cost_lid_pp=150,
                        om_lid_bc=0.08,
                        om_lid_pp=0.04,
                        discount_rate=0.02,
                        project_lifetime=30,
                    )

                # Run SWMM and get Qmax for this scenario
                subprocess.run(
                    [swmm_exe_path, inp_sc, rpt_sc],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                q_here = get_system_qmax_from_rpt(rpt_sc)
                q_s_list.append(float(q_here))

            q_mean = float(np.mean(q_s_list))
            q = q_mean
            c = float(lcc_val)
            '''
            if RAIN_TS_LIST:
                q_s_list = []
                lcc_val = None

                base_inp, ext_inp = os.path.splitext(inp_file)
                base_rpt, ext_rpt = os.path.splitext(rpt_file)

                # LCC independent of rainfall → compute once
                if lcc_val is None:
                    ts0 = RAIN_TS_LIST[0]
                    # --- make a UNIQUE inp for LCC to avoid Windows file-lock conflicts --- ** ensemble
                    pid = os.getpid()
                    uid = f"init{i}_pid{pid}_{int(time.time()*1000)}"
                    #inp_tmp = f"{base_inp}_{ts0}{ext_inp}"
                    inp_tmp = f"{base_inp}_{ts0}_LCC_{uid}{ext_inp}"

                    generate_inp_file(
                        template_inp_path,
                        inp_tmp,
                        pipe_enabled,
                        pipe_diams,
                        conduits_info,
                        lid_vars,
                        loops,
                        subcatch_info,
                        coords_info,
                        rain_gage=ts0,
                    )

                    lcc_val = compute_LCC(
                        inp_tmp, cost_per_diameter, lid_vars,
                        cost_lid_bc=200, cost_lid_pp=150,
                        om_lid_bc=0.08, om_lid_pp=0.04,
                        discount_rate=0.02, project_lifetime=30
                    )

                # === Parallel SWMM evaluation ===
                # **ensemble
                run_uid = f"init{i}_pid{os.getpid()}_{int(time.time()*1000)}"
                args_list = [
                    (
                        ts_name,
                        base_inp, ext_inp,
                        base_rpt, ext_rpt,
                        pipe_enabled,
                        pipe_diams,
                        lid_vars,
                        conduits_info,
                        loops,
                        subcatch_info,
                        coords_info,
                        template_inp_path,
                        swmm_exe_path, run_uid
                    )
                    for ts_name in RAIN_TS_LIST
                ]

                # Allocate processes based on CPU core count
                max_workers = min(len(RAIN_TS_LIST), os.cpu_count() or 4)

                with ProcessPoolExecutor(max_workers=max_workers) as ex:
                    q_s_list = list(ex.map(eval_one_scenario_parallel, args_list))

                q_mean = float(np.mean(q_s_list))
                q = q_mean
                c = float(lcc_val)
            print(
                f"[init {i}] TIMESERIES = {RAIN_TS_LIST}, "
                f"Qmax = {q_s_list}, mean={q_mean:.2f}, LCC={c:.2f}"
            )
        else:
            # fallback: original single-scenario behaviour
            generate_inp_file(
                template_inp_path,
                inp_file,
                pipe_enabled,
                pipe_diams,
                conduits_info,
                lid_vars,
                loops,
                subcatch_info,
                coords_info,
            )

            c = compute_LCC(
                inp_file,
                cost_per_diameter,
                lid_vars,
                cost_lid_bc=200,
                cost_lid_pp=150,
                om_lid_bc=0.08,
                om_lid_pp=0.04,
                discount_rate=0.02,
                project_lifetime=30,
            )

            subprocess.run(
                [swmm_exe_path, inp_file, rpt_file],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
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
        #y_init.append(w2 * qn_ini + w1 * cn_ini)
        y_init.append(w2_local * qn_ini + w1_local * cn_ini)

        qmax_ini_list.append(qn_ini)
        lcc_ini_list.append(cn_ini)
        print(
            f"Sample: qmax={q:.2f}, lcc={c:.2f}, normalized y={y_init[-1]:.4f}, "
            f"Q_MAX_MIN={Q_MAX_MIN:.4f}, Q_MAX_MAX={Q_MAX_MAX:.4f}, "
            f"LCC_MIN={LCC_MIN:.4f}, LCC_MAX={LCC_MAX:.4f}"
        )

    X_init = X_init.astype(np.float32)
    X_init = X_init[idx]
    y_init = np.array(y_init, dtype=np.float32)
    y_init = y_init[idx]
    qmax_ini_list = np.array(qmax_ini_list, dtype=np.float32)[idx]
    lcc_ini_list = np.array(lcc_ini_list, dtype=np.float32)[idx]
    print(f"[init] total: {time.time() - tsampcheck:.3f}s")
    print("=== Debug: y_init summary ===")
    print("  y_init shape:", y_init.shape)
    print("  unique values in y_init:", np.unique(y_init))
    print("  min, max, mean:", y_init.min(), y_init.max(), y_init.mean())

    # ========================
    # Bayesian Optimization Loop (Warm-start GP + y standardization)
    # ========================
    # De-duplicate initial training set
    X_init, y_init = dedupe_xy(X_init, y_init)
    y_mean = float(y_init.mean())
    y_std = float(max(y_init.std(), 1e-3))
    y_norm = (y_init - y_mean) / y_std
    train_x = torch.tensor(reorder_batch(X_init), dtype=torch.float)
    train_y = torch.tensor(y_norm, dtype=torch.float)

    rf_state = None
    if SURROGATE == "rf":
        rf_state = train_rf_regressor(
            X_init.astype(np.float32), y_init.astype(np.float32)
        )
    else:
        likelihood = gpytorch.likelihoods.GaussianLikelihood(
            noise_constraint=GreaterThan(1e-2)
        )
        likelihood.register_prior(
            "noise_prior", GammaPrior(2.0, 200.0), "noise"
        )
        model = MixedGPModel(train_x, train_y, likelihood, input_dim=D)
        with torch.no_grad():
            likelihood.noise = torch.tensor(0.01)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
        optimizer_gp = torch.optim.Adam(
            list(model.parameters()) + list(likelihood.parameters()), lr=0.01
        )

        def train_gp(train_x_t, train_y_t, max_iter=120, early_patience=25):
            model.train()
            likelihood.train()
            best = 1e9
            stale = 0
            for j in range(max_iter):
                optimizer_gp.zero_grad()
                # Added try/except to avoid NotPSD
                try:
                    with gpytorch.settings.cholesky_jitter(1e-2), gpytorch.settings.max_cholesky_size(2000):
                        out = model(train_x_t)
                        loss = -mll(out, train_y_t)
                except Exception as e:
                    # Last-resort stabilizer for NotPSD
                    if "NotPSDError" in str(type(e)) or "positive definite" in str(e):
                        with torch.no_grad():
                            likelihood.noise = torch.clamp(likelihood.noise * 2.0, min=1e-2, max=0.2)
                        print(f"[WarmGP] WARNING NotPSD -> increased noise to {likelihood.noise.item():.4f} and continue")
                        continue
                    else:
                        raise

                loss.backward()
                optimizer_gp.step()
                improved = loss.item() + 1e-4 < best
                if improved:
                    best = loss.item()
                    stale = 0
                else:
                    stale += 1
                if (j + 1) % 20 == 0 or j == 0:
                    print(
                        f"[WarmGP] iter {j+1}/{max_iter} "
                        f"loss={loss.item():.5f} noise={likelihood.noise.item():.4f}"
                    )
                if stale >= early_patience:
                    print(f"[WarmGP] early stop at {j+1}, best={best:.5f}")
                    break
            model.eval()
            likelihood.eval()

        print(">>> Initial GP hyperparameter training (warm start)")
        train_gp(train_x, train_y, max_iter=150, early_patience=30)

    if SURROGATE == "rf":
        print(">>> Initial RF training done")

    times = []  # for the time plot
    t_train_list = []
    t_cand_list = []
    t_ei_list = []
    t_feas_list = []
    t_inp_list = []
    t_swmm_list = []
    t_qmax_list = []
    t_obj_list = []
    t_iter_list = []

    lcc_list = []
    qmax_list = []
    y_next_list = []
    iter_time = time.time()

    for it in range(N_initial, N_initial + max_iter):
        t_iter_start = time.time()  # [TIMING] 单次迭代总耗时

    # b) Incremental surrogate update
        t_train_start = time.time()
        print(">>> Incremental surrogate update (", SURROGATE, ")")
        # De-duplicate and re-standardize
        X_init, y_init = dedupe_xy(X_init, y_init)
        y_mean = float(y_init.mean())
        y_std = float(max(y_init.std(), 1e-3))
        y_norm = (y_init - y_mean) / y_std

        if SURROGATE == "rf":
            rf_state = train_rf_regressor(
                X_init.astype(np.float32), y_init.astype(np.float32)
            )
        else:
            train_x = torch.tensor(reorder_batch(X_init), dtype=torch.float)
            train_y = torch.tensor(y_norm, dtype=torch.float)
            model.set_train_data(inputs=train_x, targets=train_y, strict=False)
            train_gp(train_x, train_y, max_iter=50, early_patience=12)

        t_train = time.time() - t_train_start
        t_train_list.append(t_train)
        if SURROGATE == "rf":
            print(f"[Iter {it}] RF (incremental) training: {t_train:.3f}s")
        else:
            print(
                f"[Iter {it}] GP (incremental) training: {t_train:.3f}s "
                f"noise={likelihood.noise.item():.4f}"
            )

    # c) Propose next point via acquisition
        t_cand_start = time.time()
        f_best = ((y_init.min() - y_mean) / y_std).item()
        xi = xi_schedule_piecewise(it, N_initial, max_iter)
        # Build hybrid candidate pool: random + local + top-hist
        n_random, n_local, n_best = 200, 500, 200
        best_x = X_init[np.argmin(y_init)]
        X_candidates = []
        for _ in range(n_random):
            x = np.random.rand(D)
            x = project_to_feasible_continuous(
                x,
                N_pipes=N_pipes,
                M_nodes=M,
                L=len(weight_list),
                node_list=node_list,
                conduits_info=conduits_info,
                mst_indices=sorted(mst_indices),
                system_outfall_nodes=system_outfall_node,
                catchment_areas=catchment_areas,
                max_lid_frac=0.10,
                force_level="mid",
            )
            X_candidates.append(x)
        for _ in range(n_local):
            noise = np.random.normal(loc=0.0, scale=0.05, size=D)
            x = np.clip(best_x + noise, 0.0, 1.0)
            x = project_to_feasible_continuous(
                x,
                N_pipes=N_pipes,
                M_nodes=M,
                L=len(weight_list),
                node_list=node_list,
                conduits_info=conduits_info,
                mst_indices=sorted(mst_indices),
                system_outfall_nodes=system_outfall_node,
                catchment_areas=catchment_areas,
                max_lid_frac=0.10,
                force_level="mid",
            )
            X_candidates.append(x)
        X_candidates.extend(X_init[np.argsort(y_init)[:n_best]])
        X_candidates = np.unique(np.round(X_candidates, 6), axis=0)

        if SURROGATE == "rf":
            X_feat = reorder_batch(X_candidates).astype(np.float32)
            rf = rf_state["rf"]
            preds = (
                np.array([est.predict(X_feat) for est in rf.estimators_])
                if hasattr(rf, "estimators_")
                else rf.predict(X_feat)[None, :]
            )
            mu = preds.mean(axis=0)
            std = preds.std(axis=0)  # variance proxy across trees
            beta = 1.5
            acq = mu - beta * std  # minimize
            last_x = reorder_for_kernel(X_init[-1])
            distances = np.linalg.norm(X_feat - last_x, axis=1)
            acq -= 1.0 * np.exp(-distances**2 / 1e-5)
            x_next_raw = X_candidates[int(np.argmin(acq))]
        else:
            X_cand_torch = torch.tensor(
                reorder_batch(X_candidates), dtype=torch.float
            )
            ei_vals = expected_improvement(X_cand_torch, f_best, model, xi=xi)
            last_x = reorder_for_kernel(X_init[-1])
            ei_np = ei_vals.detach().cpu().numpy().ravel()
            distances = np.linalg.norm(
                reorder_batch(X_candidates) - last_x, axis=1
            )
            ei_np -= 1.0 * np.exp(-distances**2 / 1e-5)
            x_next_raw = X_candidates[np.argmax(ei_np)]

    # Lightweight neighborhood selection (1 true eval)
        k_local = 64
        neigh = []
        for _ in range(k_local):
            z = np.clip(
                x_next_raw + np.random.normal(0, 0.02, size=x_next_raw.shape),
                0,
                1,
            )
            z = project_to_feasible_continuous(
                z,
                N_pipes=N_pipes,
                M_nodes=M,
                L=len(weight_list),
                node_list=node_list,
                conduits_info=conduits_info,
                mst_indices=sorted(mst_indices),
                system_outfall_nodes=system_outfall_node,
                catchment_areas=catchment_areas,
                max_lid_frac=0.10,
                force_level="mid",
            )
            neigh.append(z)

        neigh = np.unique(np.round(np.asarray(neigh), 6), axis=0)
        if neigh.shape[0] == 0:
            x_eval = x_next_raw.copy()
        else:
            if SURROGATE == "rf":
                Xn = reorder_batch(neigh).astype(np.float32)
                rf = rf_state["rf"]
                preds = (
                    np.array([est.predict(Xn) for est in rf.estimators_])
                    if hasattr(rf, "estimators_")
                    else rf.predict(Xn)[None, :]
                )
                mu = preds.mean(axis=0)
                std = preds.std(axis=0)
                beta = 1.5
                ucb = mu - beta * std
                x_eval = neigh[int(np.argmin(ucb))]
            else:
                with torch.no_grad():
                    Nt = torch.tensor(
                        reorder_batch(neigh), dtype=torch.float
                    )
                    post = model(Nt)
                    mu = post.mean
                    sigma = post.variance.clamp_min(1e-12).sqrt()
                    beta = 1.5
                    ucb = mu - beta * sigma
                    x_eval = neigh[int(torch.argmin(ucb))]

    # d) Evaluate x_next via SWMM
        x_next = project_to_feasible_continuous(
            x_eval,
            N_pipes=N_pipes,
            M_nodes=M,
            L=len(weight_list),
            node_list=node_list,
            conduits_info=conduits_info,
            mst_indices=mst_indices,
            system_outfall_nodes=system_outfall_node,
            catchment_areas=catchment_areas,
            max_lid_frac=0.10,
            force_level="mid",
        )
        pipe_enabled, pipe_diameters, outlet_selection, lid_vars = vector_to_design(
            x_next,
            N_pipes=N_pipes,
            M_nodes=M,
            weight_list=weight_list,
            enforce_monotone=True,
            edge_to_pipe_idx=edge_to_pipe_idx,
            fixed_outlet_mask=OUTLET_FIXED_MASK,
        )

    # enforce MST edges enabled
        for idx_mst in mst_indices:
            pipe_enabled[idx_mst] = 1
        
    # timing accumulators for this iteration
        t_inp_total = 0.0
        t_lcc_total = 0.0
        t_qmax_total = 0.0

        if RAIN_TS_LIST:
            q_s_list = []
            lcc_val = None

            base_inp, ext_inp = os.path.splitext(inp_file)
            base_rpt, ext_rpt = os.path.splitext(rpt_file)

            '''
            for ts_name in RAIN_TS_LIST:
                inp_sc = f"{base_inp}_{ts_name}{ext_inp}"
                rpt_sc = f"{base_rpt}_{ts_name}{ext_rpt}"

                # (1) Build .inp for this scenario
                t0 = time.time()
                generate_inp_file(
                    template_inp_path,
                    inp_sc,
                    pipe_enabled,
                    pipe_diameters,
                    conduits_info,
                    lid_vars,
                    loops,
                    subcatch_info,
                    coords_info,
                    rain_gage=ts_name,
                )
                debug_print_raingages(inp_sc)
                t_inp_total += time.time() - t0

                # (2) LCC (only once per design, independent of rainfall)
                if lcc_val is None:
                    t1 = time.time()
                    lcc_val = compute_LCC(
                        inp_sc,
                        cost_per_diameter,
                        lid_vars,
                        cost_lid_bc=200,
                        cost_lid_pp=150,
                        om_lid_bc=0.08,
                        om_lid_pp=0.04,
                        discount_rate=0.02,
                        project_lifetime=30,
                    )
                    t_lcc_total += time.time() - t1

                # (3) SWMM run + Qmax for this scenario
                t2 = time.time()
                subprocess.run(
                    [swmm_exe_path, inp_sc, rpt_sc],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                q_here = get_system_qmax_from_rpt(rpt_sc)
                t_qmax_total += time.time() - t2

                q_s_list.append(float(q_here))

            t_inp_list.append(t_inp_total)
            t_swmm_list.append(t_lcc_total)
            t_qmax_list.append(t_qmax_total)
            print(
                f"[Iter {it}] generate_inp_file (all models): {t_inp_total:.3f}s"
            )
            print(f"[Iter {it}] LCC compute: {t_lcc_total:.3f}s")
            print(
                f"[Iter {it}] SWMM+Qmax (all models): {t_qmax_total:.3f}s"
            )

            qmax_raw = float(np.mean(q_s_list))  # mean Qmax over RAIN_TS_LIST
            lcc = float(lcc_val)
            '''
            if RAIN_TS_LIST:
                q_s_list = []
                lcc_val = None

                base_inp, ext_inp = os.path.splitext(inp_file)
                base_rpt, ext_rpt = os.path.splitext(rpt_file)

                # LCC independent of rainfall → compute once
                if lcc_val is None:
                    ts0 = RAIN_TS_LIST[0]
                    # --- make a UNIQUE inp for LCC to avoid Windows file-lock conflicts --- ** for ensemble
                    pid = os.getpid()
                    uid = f"init{i}_pid{pid}_{int(time.time()*1000)}"
                    #inp_tmp = f"{base_inp}_{ts0}{ext_inp}"
                    inp_tmp = f"{base_inp}_{ts0}_LCC_{uid}{ext_inp}"

                    generate_inp_file(
                        template_inp_path,
                        inp_tmp,
                        pipe_enabled,
                        pipe_diameters,
                        conduits_info,
                        lid_vars,
                        loops,
                        subcatch_info,
                        coords_info,
                        rain_gage=ts0,
                    )

                    lcc_val = compute_LCC(
                        inp_tmp, cost_per_diameter, lid_vars,
                        cost_lid_bc=200, cost_lid_pp=150,
                        om_lid_bc=0.08, om_lid_pp=0.04,
                        discount_rate=0.02, project_lifetime=30
                    )

                # === Parallel SWMM run ===
                # ** for ensemble
                run_uid = f"init{i}_pid{os.getpid()}_{int(time.time()*1000)}"
                args_list = [
                    (
                        ts_name,
                        base_inp, ext_inp,
                        base_rpt, ext_rpt,
                        pipe_enabled,
                        pipe_diameters,
                        lid_vars,
                        conduits_info,
                        loops,
                        subcatch_info,
                        coords_info,
                        template_inp_path,
                        swmm_exe_path,
                        run_uid
                    )
                    for ts_name in RAIN_TS_LIST
                ]

                max_workers = min(len(RAIN_TS_LIST), os.cpu_count() or 4)

                with ProcessPoolExecutor(max_workers=max_workers) as ex:
                    q_s_list = list(ex.map(eval_one_scenario_parallel, args_list))

                qmax_raw = float(np.mean(q_s_list))
                lcc = float(lcc_val)
            print(
                f"[Iter {it}] TIMESERIES = {RAIN_TS_LIST}, "
                f"Qmax list = {q_s_list}, mean={qmax_raw:.2f}, LCC={lcc:.2f}"
            )
        else:
            # ---- Original single-scenario behaviour ----
            t_inp_start = time.time()
            generate_inp_file(
                template_inp_path,
                inp_file,
                pipe_enabled,
                pipe_diameters,
                conduits_info,
                lid_vars,
                loops,
                subcatch_info,
                coords_info,
            )
            t_inp = time.time() - t_inp_start
            t_inp_list.append(t_inp)
            print(f"[Iter {it}] generate_inp_file: {t_inp:.3f}s")

            t_swmm_start = time.time()
            lcc = compute_LCC(
                inp_file,
                cost_per_diameter,
                lid_vars,
                cost_lid_bc=200,
                cost_lid_pp=150,
                om_lid_bc=0.08,
                om_lid_pp=0.04,
                discount_rate=0.02,
                project_lifetime=30,
            )
            t_lcc = time.time() - t_swmm_start
            t_swmm_list.append(t_lcc)
            print(f"[Iter {it}] LCC compute: {t_lcc:.3f}s")

            tqmax_start = time.time()
            subprocess.run(
                [swmm_exe_path, inp_file, rpt_file],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            qmax_raw = get_system_qmax_from_rpt(rpt_file)
            t_qmax = time.time() - tqmax_start
            t_qmax_list.append(t_qmax)
            print(f"[Iter {it}] parse RPT: {t_qmax:.3f}s")

    # (4) Normalization and scalar objective (same as before)
        t_obj_start = time.time()
        qn = (qmax_raw - Q_MAX_MIN) / (Q_MAX_MAX - Q_MAX_MIN + EPS)
        cn = (lcc - LCC_MIN) / (LCC_MAX - LCC_MIN + EPS)
        #y_next = w2 * qn + w1 * cn
        y_next = w2_local * qn + w1_local * cn
        t_obj = time.time() - t_obj_start
        t_obj_list.append(t_obj)
        print(
            f"[Iter {it}] objective: {t_obj:.3f}s "
            f"(y={y_next:.4f}, qn={qn:.4f}, cn={cn:.4f})"
        )

        y_next_list.append(y_next)
        lcc_list.append(cn)
        qmax_list.append(qn)
        times.append(time.time() - iter_time)

    # e) Append new sample
        X_init = np.vstack([X_init, x_next[None, :]])
        y_init = np.append(y_init, np.float32(y_next))
        qmax_ini_list = np.append(qmax_ini_list, np.float32(qn))
        lcc_ini_list = np.append(lcc_ini_list, np.float32(cn))

        t_iter = time.time() - t_iter_start
        t_iter_list.append(t_iter)
        print(f"[Iter {it}] TOTAL iteration time: {t_iter:.3f}s")

        # Convert to a Python list for best-so-far
        y_history = list(y_init)

    # Compute best-so-far sequence
    best_so_far = [min(y_history[:i + 1]) for i in range(len(y_history))]
    # Build full best-so-far for q and c across initial + BO, then slice BO segment
    q_best_full = np.minimum.accumulate(qmax_ini_list)
    c_best_full = np.minimum.accumulate(lcc_ini_list)
    qmax_best_so_far = q_best_full[N_initial:N_initial + len(qmax_list)]
    lcc_best_so_far = c_best_full[N_initial:N_initial + len(lcc_list)]

    mean_y_init = float(np.mean(y_init[:N_initial]))
    mean_qmax_init = float(np.mean(qmax_ini_list[:N_initial]))
    mean_lcc_init = float(np.mean(lcc_ini_list[:N_initial]))

    print("Decoded outlet selection:", outlet_selection)
    print("Objective value at x_next:", y_next)

    bo_best = best_so_far[N_initial:N_initial + max_iter]

    # ==== PLOTTING ====
    if plot:
        plt.subplot(1, 3, 1)
        plt.plot(
            np.arange(1, len(y_init) + 1)[:N_initial],
            y_init[:N_initial],
            marker="d",
            linestyle="-",
            label="Initial samples (y)",
        )
        plt.plot(
            range(1, max_iter + 1),
            bo_best,
            marker="o",
            label="BO per-iteration y",
        )
        plt.axhline(
            mean_y_init,
            ls="--",
            alpha=0.6,
            color="C1",
            label=f"Init mean = {mean_y_init:.3f}",
        )
        plt.scatter(
            N_initial,
            mean_y_init,
            s=90,
            marker="X",
            color="C1",
            zorder=5,
        )
        plt.xlabel("Iteration")
        plt.ylabel("Best Objective Value So Far")
        plt.title("Bayesian Optimization Convergence Curve")
        plt.tight_layout()

    iterations = np.arange(1, len(qmax_list) + 1)

    if plot:
        plt.subplot(1, 3, 2)
        plt.plot(
            iterations,
            lcc_best_so_far,
            marker="s",
            linestyle="--",
            label="LCC",
        )
        plt.axhline(
            mean_lcc_init,
            ls="--",
            alpha=0.6,
            color="C2",
            label=f"Init mean = {mean_lcc_init:.3f}",
        )
        plt.scatter(
            N_initial,
            mean_lcc_init,
            s=90,
            marker="X",
            color="C2",
            zorder=5,
        )
        plt.xlabel("Iteration")
        plt.ylabel("Normalized LCC (0-1)")
        plt.title("LCC Over Iterations")
        plt.grid(True)
        plt.tight_layout()

    if plot:
        plt.subplot(1, 3, 3)
        plt.plot(
            iterations,
            qmax_best_so_far,
            marker="^",
            linestyle="--",
            label="Qmax",
        )
        plt.axhline(
            mean_qmax_init,
            ls="--",
            alpha=0.6,
            color="C3",
            label=f"Init mean = {mean_qmax_init:.3f}",
        )
        plt.scatter(
            N_initial,
            mean_qmax_init,
            s=90,
            marker="X",
            color="C3",
            zorder=5,
        )
        plt.xlabel("Iteration")
        plt.ylabel("Normalized qmax (0-1)")
        plt.title("qmax Over Iterations")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    if plot:
        plt.figure(figsize=(6, 4))
        plt.plot(times, bo_best, marker="o", linestyle="-")
        plt.xlabel("Elapsed Time (s)")
        plt.ylabel("Best Objective So Far")
        plt.title("Convergence vs Wall-Clock Time")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    
    try:
    # Iteration indices: 1, 2, ..., len(bo_best)
        iters = np.arange(1, len(bo_best) + 1, dtype=int)

    # Results directory: follow multi-run convention
        results_dir_env = os.environ.get("RESULTS_DIR", "results")
        if results_dir_env == "SCRIPT_DIR":
            results_dir = os.path.dirname(os.path.abspath(__file__))
        else:
            results_dir = results_dir_env
        os.makedirs(results_dir, exist_ok=True)

    # Timestamp to avoid overwriting
        ts = time.strftime("%Y%m%d_%H%M%S")

    # Note: lcc_best_so_far / qmax_best_so_far are
    # “best-so-far per generation (normalized)”, lengths aligned with bo_best / times
        df = pd.DataFrame({
            "iter": iters,
            "y_best": np.asarray(bo_best, dtype=float),
            "lcc_best": np.asarray(lcc_best_so_far, dtype=float),
            "qmax_best": np.asarray(qmax_best_so_far, dtype=float),
            "elapsed_time_s": np.asarray(times, dtype=float),
        })

    #csv_path = os.path.join(results_dir, f"bo_single_run_{ts}.csv")
        suffix = f"_{tag}" if tag else ""
        csv_path = os.path.join(results_dir, f"bo_single_run_{ts}{suffix}.csv")
        df.to_csv(csv_path, index=False)
        print(f"[Bayesian_Opt] Single-run history saved to: {csv_path}")
    except Exception as e:
    # Does not affect main flow; warn only, do not raise
        print(f"[Bayesian_Opt] WARNING: failed to save single-run CSV: {e}")

    return {
        "y_best": np.array(bo_best, dtype=float),
        "lcc_best": np.array(lcc_best_so_far, dtype=float),
        "qmax_best": np.array(qmax_best_so_far, dtype=float),
        "times": np.array(times, dtype=float),
    }

        
if __name__ == "__main__":
    template_inp_path = "./AhvazNull.inp"
    swmm_exe_path = "./swmm5.exe" 
    inp_file = "BS_opt.inp"
    rpt_file = "BS_opt.rpt"
    multi = int(os.environ.get("BO_MULTI", "1"))
    def run_weight_sweep(inp_file, rpt_file, seed=2025, plot_each=False):
        """
        Run BO 9 times with w1:w2 = 1:9 ... 9:1 under the SAME scenario-averaged evaluation.
        Save:
        (1) per-run single CSV (already inside Bayesian_Opt)
        (2) one summary CSV for all weights
        (3) one combined convergence plot (9 curves)
        """
        weight_pairs = [(i/10.0, (10-i)/10.0) for i in range(1, 10)]  # (w1, w2)
        all_runs = []  # list of dicts with iter/y/lcc/q/time + weights

        # results directory (same convention as your script)
        results_dir_env = os.environ.get("RESULTS_DIR", "results")
        if results_dir_env == "SCRIPT_DIR":
            results_dir = os.path.dirname(os.path.abspath(__file__))
        else:
            results_dir = results_dir_env
        os.makedirs(results_dir, exist_ok=True)

        ts = time.strftime("%Y%m%d_%H%M%S")

        for k, (w1v, w2v) in enumerate(weight_pairs, start=1):
            tag = f"w{int(w1v*10)}_{int(w2v*10)}"
            print("\n" + "="*70)
            print(f"[SWEEP {k}/9] Running BO with w1:w2 = {int(w1v*10)}:{int(w2v*10)}")
            print("="*70)

            # IMPORTANT: different seeds per weight to avoid identical Sobol jitter etc.
            res = Bayesian_Opt(
                inp_file, rpt_file,
                seed=seed + k,
                plot=plot_each,              # optionally show each run's 3-panel
                w1_override=w1v,
                w2_override=w2v,
                tag=tag
            )

            # Align lengths
            y_best = np.asarray(res["y_best"], dtype=float)
            l_best = np.asarray(res["lcc_best"], dtype=float)
            q_best = np.asarray(res["qmax_best"], dtype=float)
            t_best = np.asarray(res["times"], dtype=float)
            iters = np.arange(1, len(y_best) + 1, dtype=int)

            df_run = pd.DataFrame({
                "iter": iters,
                "y_best": y_best,
                "lcc_best": l_best,
                "qmax_best": q_best,
                "elapsed_time_s": t_best,
                "w1_cost": w1v,
                "w2_flow": w2v,
                "tag": tag,
            })
            all_runs.append(df_run)

        # ---- (2) save one long-format summary CSV ----
        df_all = pd.concat(all_runs, axis=0, ignore_index=True)
        csv_summary = os.path.join(results_dir, f"bo_weight_sweep_{ts}.csv")
        df_all.to_csv(csv_summary, index=False)
        print(f"[SWEEP] Summary CSV saved to: {csv_summary}")

        # ---- (3) combined convergence plot: 9 curves in one figure ----
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 4))

        for df_run in all_runs:
            tag = df_run["tag"].iloc[0]
            w1v = df_run["w1_cost"].iloc[0]
            w2v = df_run["w2_flow"].iloc[0]
            plt.plot(df_run["iter"].values, df_run["y_best"].values, label=f"{int(w1v*10)}:{int(w2v*10)}")

        plt.xlabel("Iteration")
        plt.ylabel("Best Objective So Far")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(title="w1:w2", ncol=3, fontsize=9)

        fig_path_png = os.path.join(results_dir, f"bo_weight_sweep_convergence_{ts}.png")
        plt.savefig(fig_path_png, dpi=300)
        print(f"[SWEEP] Combined convergence figure saved to: {fig_path_png}")

        plt.show()
        
        # ============================================================
        # Build WIDE-FORMAT summary CSV (paper/plot friendly)
        # ============================================================

        # Determine maximum iteration length across all runs
        max_iter = max(df["iter"].max() for df in all_runs)
        iter_index = np.arange(1, max_iter + 1, dtype=int)

        wide = pd.DataFrame({"iter": iter_index})

        for df_run in all_runs:
            tag = df_run["tag"].iloc[0]        # e.g. "w1_9"
            w1v = int(df_run["w1_cost"].iloc[0] * 10)
            w2v = int(df_run["w2_flow"].iloc[0] * 10)
            suffix = f"{w1v}_{w2v}"

            df_tmp = df_run.set_index("iter")

            wide[f"y_{suffix}"] = df_tmp["y_best"].reindex(iter_index).values
            wide[f"lcc_{suffix}"] = df_tmp["lcc_best"].reindex(iter_index).values
            wide[f"qmax_{suffix}"] = df_tmp["qmax_best"].reindex(iter_index).values
            wide[f"time_{suffix}"] = df_tmp["elapsed_time_s"].reindex(iter_index).values

        csv_wide = os.path.join(results_dir, f"bo_weight_sweep_wide_{ts}.csv")
        wide.to_csv(csv_wide, index=False)

        print(f"[SWEEP] Wide-format CSV saved to: {csv_wide}")


        return csv_summary, fig_path_png

    if multi > 1:
        gens = max_iter
        base_seed = 2025
        y_runs = np.zeros((multi, gens))
        l_runs = np.zeros((multi, gens))
        q_runs = np.zeros((multi, gens))
        t_runs = np.zeros((multi, gens))
        for r in range(multi):
            res = Bayesian_Opt(inp_file, rpt_file, seed=base_seed + r, plot=False)
            y_runs[r,:] = res['y_best'][:gens]
            l_runs[r,:] = res['lcc_best'][:gens]
            q_runs[r,:] = res['qmax_best'][:gens]
            t_runs[r, :gens] = res['times'][:gens]
        # Results directory handling
        results_dir_env = os.environ.get("RESULTS_DIR", "results")
        if results_dir_env == "SCRIPT_DIR":
            results_dir = os.path.dirname(os.path.abspath(__file__))
        else:
            results_dir = results_dir_env
        os.makedirs(results_dir, exist_ok=True)

        # Save per-run and stats CSVs
        ts = time.strftime("%Y%m%d_%H%M%S")
        x = np.arange(1, gens+1)

        def save_runs_stats(prefix: str, runs: np.ndarray):
            runs = np.asarray(runs, dtype=float)
            # Per-run CSV: gen + each run column
            runs_df = pd.DataFrame(
                np.column_stack([x, runs.T]),
                columns=["gen"] + [f"run_{i+1}" for i in range(runs.shape[0])]
            )
            runs_df.to_csv(os.path.join(results_dir, f"bo_{prefix}_runs_{ts}.csv"), index=False)
            # Stats CSV
            m = runs.mean(axis=0); s = runs.std(axis=0)
            mn = runs.min(axis=0); mx = runs.max(axis=0)
            stats_df = pd.DataFrame({
                "gen": x,
                "mean": m,
                "std": s,
                "min": mn,
                "max": mx,
            })
            stats_df.to_csv(os.path.join(results_dir, f"bo_{prefix}_stats_{ts}.csv"), index=False)

        save_runs_stats("y", y_runs)
        save_runs_stats("l", l_runs)
        save_runs_stats("q", q_runs)
        save_runs_stats("time", t_runs)

        # Plot with mean ±1σ shaded, plus min/max overlays
        plt.figure(figsize=(12,4))

        # Objective subplot
        plt.subplot(1,3,1)
        y_m = y_runs.mean(axis=0); y_s = y_runs.std(axis=0)
        y_lo = np.minimum(y_m - y_s, y_m + y_s)
        y_hi = np.maximum(y_m - y_s, y_m + y_s)
        # Set axis then clip fill to range
        plt.ylim(0, 0.6)
        ax_ymin, ax_ymax = 0.0, 0.6
        plt.fill_between(x, np.clip(y_lo, ax_ymin, ax_ymax), np.clip(y_hi, ax_ymin, ax_ymax),
                         color='C0', alpha=0.35, label='±1σ', zorder=1)
        plt.plot(x, np.clip(y_m, ax_ymin, ax_ymax), color='C0', label='Objective mean', zorder=2)
        plt.plot(x, np.clip(y_runs.min(axis=0), ax_ymin, ax_ymax), color='C0', linestyle='--', alpha=0.9, label='min', zorder=3)
        plt.plot(x, np.clip(y_runs.max(axis=0), ax_ymin, ax_ymax), color='C0', linestyle='--', alpha=0.9, label='max', zorder=3)
        plt.title('Bayesian Optimization Convergence Curve')
        plt.xlabel('Iteration')
        plt.ylabel('Best Objective Value So Far')
        plt.grid(True); plt.legend()

        # LCC subplot
        plt.subplot(1,3,2)
        l_m = l_runs.mean(axis=0); l_s = l_runs.std(axis=0)
        l_lo = np.minimum(l_m - l_s, l_m + l_s)
        l_hi = np.maximum(l_m - l_s, l_m + l_s)
        plt.ylim(0, 1.0)
        ax_ymin, ax_ymax = 0.0, 1.0
        plt.fill_between(x, np.clip(l_lo, ax_ymin, ax_ymax), np.clip(l_hi, ax_ymin, ax_ymax),
                         color='green', alpha=0.35, label='±1σ', zorder=1)
        plt.plot(x, np.clip(l_m, ax_ymin, ax_ymax), color='green', label='LCC mean', zorder=2)
        plt.plot(x, np.clip(l_runs.min(axis=0), ax_ymin, ax_ymax), color='green', linestyle='--', alpha=0.9, label='min', zorder=3)
        plt.plot(x, np.clip(l_runs.max(axis=0), ax_ymin, ax_ymax), color='green', linestyle='--', alpha=0.9, label='max', zorder=3)
        plt.title('Avg LCC (norm)')
        plt.xlabel('Iteration'); plt.ylabel('LCC (norm)')
        plt.grid(True); plt.legend()

        # Qmax subplot
        plt.subplot(1,3,3)
        q_m = q_runs.mean(axis=0); q_s = q_runs.std(axis=0)
        q_lo = np.minimum(q_m - q_s, q_m + q_s)
        q_hi = np.maximum(q_m - q_s, q_m + q_s)
        plt.ylim(0, 0.6)
        ax_ymin, ax_ymax = 0.0, 0.6
        plt.fill_between(x, np.clip(q_lo, ax_ymin, ax_ymax), np.clip(q_hi, ax_ymin, ax_ymax),
                         color='orange', alpha=0.35, label='±1σ', zorder=1)
        plt.plot(x, np.clip(q_m, ax_ymin, ax_ymax), color='orange', label='Qmax mean', zorder=2)
        plt.plot(x, np.clip(q_runs.min(axis=0), ax_ymin, ax_ymax), color='orange', linestyle='--', alpha=0.9, label='min', zorder=3)
        plt.plot(x, np.clip(q_runs.max(axis=0), ax_ymin, ax_ymax), color='orange', linestyle='--', alpha=0.9, label='max', zorder=3)
        plt.title('Avg Qmax (norm)')
        plt.xlabel('Iteration'); plt.ylabel('Qmax (norm)')
        plt.grid(True); plt.legend()

        plt.tight_layout(); plt.show()
    else:
        #Bayesian_Opt(inp_file, rpt_file, plot=True)
        # weight sweep 1:9 ... 9:1
        run_weight_sweep(inp_file, rpt_file, seed=2025, plot_each=False)



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