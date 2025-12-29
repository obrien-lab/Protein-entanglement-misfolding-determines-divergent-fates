#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects
import pyemma as pem
import parmed as pmd
import mdtraj as mdt
import msmtools
import networkx as nx

matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'medium'
matplotlib.rcParams['ytick.labelsize'] = 'medium'
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['figure.dpi'] = 600

################################# Arguments ###################################
n_traj = 50
num_rep = 1
mutant_type_list = ['A6NDU8']

co_dir = '/storage/home/yuj179/mygroup/protein_Ubq/continuous_synthesis/A6NDU8/1-50/analysis/'
co_msm_data_file = co_dir+'msm_data.npz'
po_dir = '/storage/home/yuj179/mygroup/protein_Ubq/post_translation/A6NDU8/1-50/analysis/'
post_msm_data_file = po_dir+'msm_data.npz'

A = [[1]]
B = ['end']

filter_threshold = 1.0
table_filter_threshold = 15

trap_mask = []
native_mask = ['P3']

color_list = np.array([[255, 127, 0],
                       [0,   0,   255]])
color_list = color_list / 255

################################# Functions ###################################
def read_pathways(path_file):
    f = open(path_file)
    lines = f.readlines()
    f.close()
    pathways = {}
    numbers = [str(i) for i in range(10)]
    for line in lines:
        line = line.strip()
        if line == '':
            continue
        elif line.endswith('pathways:'):
            mutant = line.split()[0]
            pathways[mutant] = []
        elif line[0] in numbers:
            words = line.split()
            pathways[mutant].append([words[2:len(words):2], float(words[0])/100])
    return pathways
    
def pathway_filter(pathways, A, B, threshold):
    sub_pathways = []
    tot_fraction = 0
    for path in pathways:
        path_0 = path[0]
        if int(path_0[0]) in A and int(path_0[-1]) in B:
            sub_pathways.append(path)
            tot_fraction += path[1]
    for i in range(len(sub_pathways)):
        sub_pathways[i][1] /= tot_fraction
    sub_pathways_sorted = sorted(sub_pathways, key=lambda x: x[1], reverse=True)
    pathways_filtered = []
    if threshold > 1:
        i = 0
        for path in sub_pathways_sorted:
            pathways_filtered.append(path)
            i += 1
            if i >= threshold:
                break
    else:
        fraction = 0
        for path in sub_pathways_sorted:
            pathways_filtered.append(path)
            fraction += path[1]
            if fraction >= threshold:
                break
    return pathways_filtered
    
def build_network(count_matrix, states_on_pathway, node_label_list):
    n_state = len(count_matrix)
    weight = np.zeros((n_state,n_state))
    for i in range(n_state):
        for j in range(n_state):
            if i != j and count_matrix[i, j] > 0:
                weight[i,j] = 1
    
    G = nx.DiGraph()
    nodes = []
    for i in range(0,n_state):
        if i not in states_on_pathway:
            continue
        nodes.append((i, {'label':node_label_list[i],
                          'size':50,}))
    G.add_nodes_from(nodes)

    edges = []
    for i in range(0,n_state):
        for j in range(0,n_state):
            if i == j or weight[i,j] == 0:
                continue
            edges.append((i, j, {'weight':weight[i,j],
                                 'viz':{'thickness':weight[i,j],
                                        'shape':'solid',
                                        'color':{'a':1,
                                                 'r':0,
                                                 'g':0,
                                                 'b':0}}}))
    G.add_edges_from(edges)
    
    return G

def get_graph_layout(G):
    n_state = len(list(G.nodes()))
    initial_pos = {}
    for idx, (node, data) in enumerate(G.nodes(data=True)):
        x = idx/(n_state-1)
        #if idx%2 == 0:
            #y = 0.5
        #else:
            #y = -0.5
        y = np.random.rand()
        initial_pos[node] = [x,y]
    print(initial_pos)
    pos = nx.kamada_kawai_layout(G, pos=initial_pos, dim=2)
    print(pos)
    fixed_node = [0]
    # find close contacts
    cutoff = 0.1
    for u in G.nodes():
        if u == 0:
            continue
        tag_fix = True
        for v in G.nodes():
            if u == v:
                continue
            distance = ((pos[u][0] - pos[v][0])**2 + (pos[u][1] - pos[v][1])**2)**0.5
            if distance <= cutoff:
                tag_fix = False
                break
        if tag_fix:
            fixed_node.append(u)
    #print(fixed_node)
    fixed_node = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 20, 21, 22]
    pos = nx.spring_layout(G, k=0.02, pos=pos, fixed=fixed_node, seed=1)
    return pos
    
def read_layout(file_name):
    G = nx.read_gml(file_name, label=None)
    pos_0 = []
    for idx, (node, data) in enumerate(G.nodes(data=True)):
        pos_0.append([data['graphics']['x'], data['graphics']['y']])
    pos_0 = np.array(pos_0)
    pos_0 = nx.rescale_layout(pos_0)
    pos = {}
    for idx, (node, data) in enumerate(G.nodes(data=True)):
        pos[node] = list(pos_0[idx])
    return(pos)

def read_layout_gexf(file_name):
    import xml.etree.ElementTree as ET
    import re
    
    nsmap = {}
    for event, (prefix, uri) in ET.iterparse(file_name, events=['start-ns']):
        # later declarations override earlier ones for same prefix
        nsmap[prefix or ''] = uri
    
    tree = ET.parse(file_name)
    root = tree.getroot()
    m = re.match(r"\{(.+)\}", root.tag)
    base_ns = m.group(1) if m else nsmap.get('', '')
    viz_ns = nsmap.get('viz')
    
    ns = {"g": base_ns, "viz": viz_ns}
    
    pos_0 = []
    for node_elem in root.findall(".//g:node", ns):
        nid = node_elem.get("id")
        viz_pos = node_elem.find("viz:position", ns)
        if viz_pos is not None:
            pos_0.append([float(viz_pos.get('x')), float(viz_pos.get('y'))])

    pos_0 = np.array(pos_0)
    pos_0 = nx.rescale_layout(pos_0)
    pos = {}
    
    for idx, node_elem in enumerate(root.findall(".//g:node", ns)):
        nid = int(node_elem.get("id"))
        pos[nid] = list(pos_0[idx])
    return(pos)


def draw_pathways(ax, G, pos):
    global image_list
    n_state = len(list(G.nodes()))
    
    for (node, data) in G.nodes(data=True):
        if data['label'] in trap_mask:
            edgecolor = 'r'
        elif data['label'] in native_mask:
            edgecolor = np.array([46, 139, 87])/255
        else:
            edgecolor = 'k'
        nx.draw_networkx_nodes(G, pos, nodelist=[node], ax=ax,
                                node_size=800, 
                                alpha=1,
                                linewidths=1.0,
                                edgecolors=edgecolor,
                                node_color='w',
                                node_shape='s')
    
    for (u,v,d) in G.edges(data=True):
        color = np.array([[d['viz']['color']['r']/255, d['viz']['color']['g']/255, d['viz']['color']['b']/255]])
        collection = nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], ax=ax,
                                            width=1*d['viz']['thickness'],
                                            edge_color=color,
                                            alpha=1,
                                            arrowsize=10,
                                            connectionstyle='arc3,rad=0.2',
                                            node_size=800)
        
    node_labels = {}
    for (node, data) in G.nodes(data=True):
        if data['label'] in trap_mask:
            node_labels[node] = data['label']+'/Trapped'
        elif data['label'] in native_mask:
            node_labels[node] = data['label']+'/Folded'
        else:
            node_labels[node] = data['label']
    label_pos = {}
    for p in pos.keys():
        label_pos[p] = list(np.array(pos[p])+np.array([0.05, 0.05]))
    label_dict = nx.draw_networkx_labels(G, label_pos, labels=node_labels, ax=ax, font_size=8,
                                         horizontalalignment='left', verticalalignment='center')
    for i in label_dict.keys():
        label_dict[i].set_path_effects([matplotlib.patheffects.withStroke(linewidth=2, foreground='w')])
        label_dict[i].set_zorder(10)
    ax.axis('off')
    
    # add annotation images
    for node in G.nodes():
        arr_img = image_list[node]
        zoom = 25/max(arr_img.shape[:2])
        
        imagebox = matplotlib.offsetbox.OffsetImage(arr_img, zoom=zoom)
        imagebox.image.axes = ax

        ab = matplotlib.offsetbox.AnnotationBbox(imagebox, pos[node],
                                                 xycoords='data', frameon=False)

        ax.add_artist(ab)
    ax.set_xlim(np.array(ax.get_xlim()) + np.abs(ax.get_xlim())*np.array([-0.05, 0.3]))
    ax.set_ylim(np.array(ax.get_ylim()) + np.abs(ax.get_ylim())*np.array([-0.2, 0.3]))
    # print(ax.get_xlim(), ax.get_ylim())
    
def show_pathways_table(ax, pathways, node_label_list):
    global table_filter_threshold
    cell_content = []
    for path in pathways:
        words = []
        for p in path[0]:
            words.append(node_label_list[int(p)-1])
        cell_content.append(['%.2f %%'%(path[1] * 100), u' \u2192 '.join(words)])
    if table_filter_threshold <= 1:
        column_text = ['Percentage', 'Top %d %% paths'%(table_filter_threshold*100)]
    else:
        column_text = ['Percentage', 'Top %d paths'%(table_filter_threshold)]
    tab = ax.table(cellText=cell_content, cellLoc='center', colLabels=column_text,
                   rowLoc='center', colLoc='center', loc='upper center')
    tab.auto_set_font_size(False)
    tab.set_fontsize(5)
    tab.auto_set_column_width(col=list(range(len(column_text))))
    cellDict = tab.get_celld()
    for cell in cellDict.keys():
        cellDict[cell].set_height(0.05)
    
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
def show_statistic_table(ax, pathways, node_label_list):
    from scipy.stats import bootstrap
    cell_content = []
    stat = {}
    for path in pathways:
        C_last = ''
        P_last = ''
        for p in path[0]:
            name = node_label_list[int(p)-1]
            if name.startswith('C'):
                C_last = name
            else:
                P_last = name
        label = u'* \u2192 %s \u2192 * \u2192 %s'%(C_last, P_last)
        if label not in stat.keys():
            stat[label] = path[1]
        else:
            stat[label] += path[1]
    sort_stat = sorted(stat.items(), key=lambda x: x[1], reverse=True)
    
    for s in sort_stat:
        data = np.zeros(n_traj)
        data[:int(n_traj*s[1])] = 1
        res = bootstrap((data,), np.mean, method='percentile')
        ci_l, ci_u = res.confidence_interval
        cell_content.append([s[0], '%.2f %% [%.2f %%, %.2f %%]'%(s[1]*100, ci_l*100, ci_u*100)])
    column_text = ['Path patterns', 'Percentage']
    tab = ax.table(cellText=cell_content, cellLoc='center', colLabels=column_text,
                   rowLoc='center', colLoc='center', loc='upper center')
    tab.auto_set_font_size(False)
    tab.set_fontsize(5)
    tab.auto_set_column_width(col=list(range(len(column_text))))
    cellDict = tab.get_celld()
    for cell in cellDict.keys():
        cellDict[cell].set_height(0.05)
    
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
def add_inset_axes(ax, x, y):
    ins = ax.inset_axes([0.5,0.4,0.3,0.3])
    ins.plot(x, y, '-k')
    ins.plot([x[0], x[-1]], [0, 0], '-r')
    return ins
    
################################## MAIN #######################################
# Read co-trans metastable dtrajs
co_meta_dtrajs = np.load(co_msm_data_file, allow_pickle=True)['meta_dtrajs']
co_n_states = 0
for i in range(len(co_meta_dtrajs)):
    if np.max(co_meta_dtrajs[i]) > co_n_states:
        co_n_states = np.max(co_meta_dtrajs[i])
co_n_states += 1
# Read post-trans metastable dtrajs
post_meta_dtrajs = np.load(post_msm_data_file, allow_pickle=True)['meta_dtrajs']
post_n_states = 0
for i in range(len(post_meta_dtrajs)):
    if np.max(post_meta_dtrajs[i]) > post_n_states:
        post_n_states = np.max(post_meta_dtrajs[i])
post_n_states += 1
node_label_list = []
image_list = []
for i in range(co_n_states):
    node_label_list.append('C%d'%(i+1))
    os.system('convert %s/state_struct/state_%d.tga ./%d.png'%(co_dir, i+1, i+1))
    img = plt.imread('%d.png'%(i+1))
    image_list.append(img)
    os.system('rm -f %d.png'%(i+1))
for i in range(post_n_states):
    node_label_list.append('P%d'%(i+1))
    os.system('convert -trim %s/state_struct/state_%d.tga ./%d.png'%(po_dir, i+1, i+1))
    img = plt.imread('%d.png'%(i+1))
    image_list.append(img)
    os.system('rm -f %d.png'%(i+1))
# Read pathways
pathways = read_pathways('./pathways.dat')
end_state_id = [int(p[0][-1]) for mutant_type in mutant_type_list for p in pathways[mutant_type]]
for i in range(len(B)):
    if B[i] == 'end':
        B[i] = end_state_id

# Draw pathways
panel_idx = list(map(chr, range(ord('a'), ord('z')+1)))
fig_width = 5
num_rows = len(A) * len(B)
row_height = fig_width*0.3
fig_height = num_rows * row_height
fig = plt.figure(figsize=(fig_width, fig_height))
gs = fig.add_gridspec(nrows=num_rows*10, ncols=2, 
                      top=0.9, bottom=0.1, left=0.1, right=0.9,
                      hspace=0.05, wspace=0.1)
i_ax = 0
for start_idx_list in A:
    for end_idx_list in B:
        pathways_filtered_combined = []
        for mutant_type in mutant_type_list:
            pathways_filtered_combined += pathway_filter(pathways[mutant_type], start_idx_list, end_idx_list, filter_threshold)
        dtrajs = []
        for path in pathways_filtered_combined:
            p = np.array(path[0], dtype=int)-1
            dtrajs.append(p)
        c_matrix = msmtools.estimation.count_matrix(dtrajs, 1).toarray()
        states_on_pathway = []
        for i in range(c_matrix.shape[0]):
            if np.sum(c_matrix[i, :]) > 0 or np.sum(c_matrix[:, i]) > 0:
                states_on_pathway.append(i)
        G = build_network(c_matrix, states_on_pathway, node_label_list)
        nx.write_gexf(G, "network_%.2f.gexf"%filter_threshold)
        #pos = get_graph_layout(G)
        pos = read_layout_gexf('network_%.2f_layout.gexf'%filter_threshold)
        G_mutant_list = []
        for mutant_type in mutant_type_list:
            pathways_filtered = pathway_filter(pathways[mutant_type], start_idx_list, end_idx_list, filter_threshold)
            dtrajs = []
            for path in pathways_filtered:
                p = np.array(path[0], dtype=int)-1
                dtrajs.append(p)
            c_matrix = msmtools.estimation.count_matrix(dtrajs, 1).toarray()
            states_on_pathway = []
            for i in range(c_matrix.shape[0]):
                if np.sum(c_matrix[i, :]) > 0 or np.sum(c_matrix[:, i]) > 0:
                    states_on_pathway.append(i)
            G_0 = build_network(c_matrix, states_on_pathway, node_label_list)
            G_mutant_list.append(G_0)
        for (u,v,d) in G.edges(data=True):
            if len(G_mutant_list) > 1:
                if (u,v) in G_mutant_list[0].edges() and (u,v) not in G_mutant_list[1].edges():
                    d['viz']['color']['r'] = color_list[0,0]*255
                    d['viz']['color']['g'] = color_list[0,1]*255
                    d['viz']['color']['b'] = color_list[0,2]*255
                elif (u,v) not in G_mutant_list[0].edges() and (u,v) in G_mutant_list[1].edges():
                    d['viz']['color']['r'] = color_list[1,0]*255
                    d['viz']['color']['g'] = color_list[1,1]*255
                    d['viz']['color']['b'] = color_list[1,2]*255

        ax = fig.add_subplot(gs[i_ax*10:(i_ax*10+10), :])
        bbox = ax.get_position()
        (x, y, width, height) = bbox.bounds
        draw_pathways(ax, G, pos)
        
        i_ax += 1
fig.savefig('pathways_%.2f.svg'%(filter_threshold))

# Draw statistic tables
table_fig = plt.figure(figsize=(4,6))
for i_ax, mutant_type in enumerate(mutant_type_list):
    ax = table_fig.add_subplot(2,1,i_ax+1)
    show_statistic_table(ax, pathways[mutant_type], node_label_list)
table_fig.savefig('pathway_statistic.svg')
