# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 12:01:09 2023

@author: u03132tk
"""
import pandas as pd
import time
from Bio import Entrez, SeqIO
from http.client import IncompleteRead
from ncbixml_blast_funcs import makeblastdb_subprocess, blastP_subprocess
from Bio.Blast import NCBIXML
import json
import os
import plotly.graph_objects as go
import networkx as nx

#parse csv to get the scaffold accessions of each hit and the start/end of the hits
def plot_graph(G, path):
    G = G.to_undirected()
    pos = nx.drawing.layout.spring_layout(G)
    nx.set_node_attributes(G, pos, 'pos')
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text
    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='<br>Network graph made with Python',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                        text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.write_html(path)
def write_fasta(output_filepath, fasta_tuple):
    with open (output_filepath, 'w') as outfile:
        for index, (defline, seq) in enumerate(fasta_tuple):
            assert defline[0] == '>'
            outfile.write(defline + '\n')
            outfile.write(seq)
            if index < len(fasta_tuple) - 1:
                outfile.write('\n')

anchor_node = 'NZ_CP042324.1'
write_raw_files = True
neighborhood_size = 40000
file_flag = '_binary'
#add option to just read in raw files rather than scraping - can rebild fasta etc from that
filter_duplicate_taxids = False #DO NOT DO THIS IF POSSIBLE  - two lividans strains had very nice bgc delineation for coelicolor bgc 1
strict_span = True
similarity_filter = 0.7#0.8
#user_query_len = 7 #or none or int. if you have any hits in query that arent required then you need to specifiy minimal number 
folders = [#(r'D:/cluster_sieve', 24),
            #/sco/co-transcribed_upper', 8),#,#mibig
           #(r'D:/ModuleMapper_ClusterSieve/sco/mibig_cotranscribed', 5)
           (r'D:/ClusterSieve', 9)
           ]#'D:/CAPTURE/hyalina/bif47'#'D:/CAPTURE/herbida/r20'##'D:/CAPTURE/herbida/r26'
email = input('please enter your email for entrez searching')
missing = []
path_index_map = {}
for folder, query_len in folders:
    
    
    protein_lengths = {}
    short = 0
    #duplicate_taxids = 0
    size_map = {}
    
    parse_file = [i for i in os.listdir(folder) if file_flag in i]
    assert len(parse_file) == 1
    parse_file = parse_file[0]
    location_data = pd.read_csv(f'{folder}/{parse_file}', sep = ',')
    #if user_query_len == None:
    #    query_len = len(location_data.columns) - 5 #organism, scaffold, start, end, score
    #else:
    #    query_len = user_query_len
    location_map = {scaffold: (start, stop) for scaffold, start, stop in zip (location_data['Scaffold'],
                                                                              location_data['Start'],
                                                                              location_data['End'])}
    
    
    #scaffold_cluster_data = pd.read_csv('D:/D637S612Q784L50/results/D637S612Q784L50_summary.txt', names = range(9))[0]
    #cluster_indexes = [i for i, v in enumerate(scaffold_cluster_data) if 'Cluster' in v]
    #scaffolds = [scaffold_cluster_data[i-2] for i in cluster_indexes]
    
    #find the taxid of every scaffold, and only keep one scaffold for each taxid (redundant sequences can be misleading)
    taxid_map = {}
    #now get the regions around each hit from a unique taxid, for analysis using clinker or blast
    accessions = location_map.keys()
    Entrez.email = email
    fasta = []
    found_taxids = set()
    cut_records = {}
    print (f'Extracting {len(accessions)} gbks...')
    for index, accession in enumerate(accessions):
        print (f'accession {accession} #{index} of {len(accessions)-1}')
        if filter_duplicate_taxids:
            #check taxid last in case one record has better size
            found = False
            for term in [accession, 'NZ_'+accession]:
                if found:
                    break
                for attempt in range(5):
                    db="nucleotide"
                    handle = Entrez.esearch(db=db, term = term)
                    taxid_record = Entrez.read(handle)
                    try:
                        gi = taxid_record["IdList"][0]
                    except IndexError:
                        continue
                    found = True
                    break
            if not found:
                missing += [(accession, f'accession: {accession}')]
                print (f'CANNOT FIND - {accession}')
                continue
            found = False
            for attempt in range(5):
                handle = Entrez.esummary(db=db, id=gi, retmode="json")
                try:
                    result = json.load(handle)["result"]
                    
                except ValueError:
                    #JSONDecodeError
                    print (f'JSON issue #{attempt} - retrying...')
                    continue
                if result['uids'] == []:
                    print ('no UIDs - retrying...')
                    continue
                found = True
                break
            if not found:
                missing += [(accession, f'gi: {gi}')]
                print (f'CANNOT FIND - {gi}')
                continue
            taxid = str(result[gi]["taxid"])
            if taxid in found_taxids:
                if accession != anchor_node:
                    print ('found taxid')
                    continue
            
                
        for attempt in range(5):
            print ('scraping')
            handle = Entrez.efetch(db = "nucleotide", 
                                    id = accession, 
                                    rettype = "gbwithparts ", 
                                    retmode = "text")   
            print ('reading')
            try:
                record = SeqIO.read(handle, 
                                    "genbank")
            except IncompleteRead:
                print ('IncompleteRead - trying again')
                continue
            found = True
            break
        if not found:
            missing += [(accession, 'incomplete read')]
        motif_start, motif_stop = min(location_map[accession]), max(location_map[accession])
        #extract motif record
        #get feature names 
        #check feature name in big record so you can get indexes to ignore - {acc : ignore1, ignore2 ...}
        motif_span = motif_stop-motif_start
        if motif_span > neighborhood_size:
            raise ValueError(rf'{motif_span} larger than neighborhood size {neighborhood_size}')
        extension = (neighborhood_size - motif_span)/2
        start = motif_start - extension#extension
        stop = motif_stop + extension#extension
        
        if stop > len(record):
            if strict_span:
                short +=1
                print (f'stop ({stop}) > length of record ({len(record)}) - skipping')
                continue
            else:
                stop = len(record)
        if start < 0:
            if strict_span:
                print (f'start ({start}) < 0 - skipping')
                short += 1
                continue
            else:
                start = 0
        #you check if a taxid has been foubnd early on - but you only add a taxid to found taxids if it is written
        if filter_duplicate_taxids:
            found_taxids.update([taxid])
            #duplicate_taxids += 1
        print (f'neighborhood = {start} -> {stop}')
        #remove bias eg from contig edges
        cut_record = record[int(start) : int(stop)]
        cut_record.annotations = record.annotations
        cut_records[record.id] = cut_record
        path = f'{folder}/{record.annotations["organism"]}_raw_{index}.gbk'
        if write_raw_files:
            #note - this will overwrite files that have same organism id.  these are only a ercord of raw files so low priority
            with open(path, 'w') as handle:
                SeqIO.write (cut_record, 
                             handle, 
                             "genbank")
            print(f'{path} written: {os.path.isfile(path)}')
        #write to file for clinker #
        # with open(f'{folder}\{record.annotations["organism"]}.gbk', 'w') as handle:
        #     SeqIO.write (cut_record, 
        #                  handle, 
        #                  "genbank")
        cds_count = 0
        size_map[accession] = 0
        for feature in cut_record.features:
            if feature.type == 'CDS':
                #if feature id not in motif span record names
                #no - you need actual names in case the homologs have drifyted 
                try:
                    fasta += [(f'>{accession}__{cds_count}', 
                                ''.join(feature.qualifiers['translation']))]
                    protein_lengths[f'{accession}_{cds_count}'] = len(''.join(feature.qualifiers['translation']))
                    size_map[accession] += 1
                    cds_count += 1

                except KeyError:
                    assert 'pseudo' in feature.qualifiers.keys()
                
        
    input_db_fasta_path = f'{folder}/db.txt' 
    db_outpath = f'{folder}/db' 
    query_user = f'{folder}/query.txt'   
    write_fasta(input_db_fasta_path, fasta)
    write_fasta(query_user, fasta)
    print (f'all v all blast p - {len(fasta)} seqs')
    makeblastdb_exe_path = 'C:/Users/u03132tk/.spyder-py3/ModuleMapper/Backend/Executables/NCBI/blast-2.10.1+/bin/makeblastdb.exe'
    # query_user = 'C:/Users/u03132tk/.spyder-py3/boundary_dataset/NZ_CP042324.1.region011.fasta'
    blastp_exe_path = 'C:/Users/u03132tk/.spyder-py3/ModuleMapper/Backend/Executables/NCBI/blast-2.10.1+/bin/blastp.exe'
    results_out_path =  f'{folder}/results.xml'  
    makeblastdb_subprocess(makeblastdb_exe_path, 
                            input_db_fasta_path, 
                            db_outpath)
    print ('made db')
    start_blastp=time.time()
    blastP_subprocess (10**-5, 
                        query_user, 
                        blastp_exe_path, 
                        results_out_path, 
                        db_outpath, 
                        4)   
    raw_results = {}#
        
    with open(results_out_path, 'r') as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        for record in blast_records:
            query_scaffold, query_index = record.query.split('__')
            if query_scaffold not in raw_results.keys():
                raw_results[query_scaffold] = {}
            if query_index not in raw_results[query_scaffold].keys():
                raw_results[query_scaffold][query_index] = {}
            for hit in record.alignments:
                hit_scaffold, hit_index = hit.hit_def.split('__')
                best_hsp = None
                for hsp in hit.hsps:
                    if hsp.align_length/protein_lengths[f'{query_scaffold}_{query_index}'] < 0.5:
                        continue
                    if best_hsp == None:
                        best_hsp = hsp
                    else:
                        if hsp.score > best_hsp.score:
                            best_hsp = hsp
                
                if hit_scaffold not in raw_results[query_scaffold][query_index].keys():
                    raw_results[query_scaffold][query_index][hit_scaffold] = {}
                raw_results[query_scaffold][query_index][hit_scaffold][hit_index] = best_hsp
    best_hits = {}
    for query_scaffold, query_proteins in raw_results.items():
        best_hits[query_scaffold] = {}
        for query_index, hit_proteins in query_proteins.items():
            best_hits[query_scaffold][query_index] = {}
            for hit_scaffold, hit_proteins in hit_proteins.items():
                best_hit = None
                best_score = 0
                for hit_index, hit_hsp in hit_proteins.items():
                    if hit_hsp == None:
                        continue
                    if hit_hsp.score > best_score:
                        best_hit = hit_index
                        best_score= hit_hsp.score
                if best_hit == None:
                    #no hsps that meet coverage lims
                    continue
                best_hits[query_scaffold][query_index][hit_scaffold] = best_hit
    
    
    
        
    reciprocal_best_hits = {}
    for query_scaffold, query_proteins in best_hits.items():
        for query_index, best_hit_proteins in query_proteins.items():
            for hit_scaffold, best_hit_protein in best_hit_proteins.items():
                if query_scaffold not in best_hits[hit_scaffold][best_hit_protein].keys():
                    continue # one direction hit
                if best_hits[hit_scaffold][best_hit_protein][query_scaffold] == query_index:
                    if query_scaffold not in reciprocal_best_hits.keys():
                        reciprocal_best_hits[query_scaffold] = {}
                    if query_index not in reciprocal_best_hits[query_scaffold].keys():
                        reciprocal_best_hits[query_scaffold][query_index] = {}
                    reciprocal_best_hits[query_scaffold][query_index][hit_scaffold] = best_hit_protein
    
    #clinker_accessions = [size_map[0][0]]
    #sort by size
    
    similarity_map = {}
    for record_accession in reciprocal_best_hits.keys():
        similarity_map[record_accession] = {}
        for other_record_accession in reciprocal_best_hits.keys():
            number_of_rbh = 0
            smallest_record_protein_count = min(size_map[record_accession], size_map[other_record_accession])
            
            for query_protein, rbh_scaffolds in reciprocal_best_hits[other_record_accession].items():
                if record_accession in rbh_scaffolds.keys():
                    number_of_rbh += 1
            similarity_map[record_accession][other_record_accession] = (number_of_rbh - query_len)/(smallest_record_protein_count - query_len)
            #if sorted(['NZ_JAJLQY010000002.1', 'NZ_RJKE01000001.1']) == sorted(record_accession, other_record_accession):
            #    print (record_accession, other_record_accession, number_of_rbh)
    
    
    ##
    edge_map = {}
    for record_accession in reciprocal_best_hits.keys():
        edge_map[record_accession] = {}
        for other_record_accession in reciprocal_best_hits.keys():
            if record_accession != other_record_accession:
                number_of_rbh = 0
                smallest_record_protein_count = min(size_map[record_accession], size_map[record_accession])
                for query_protein, rbh_scaffolds in reciprocal_best_hits[other_record_accession].items():
                    if record_accession in rbh_scaffolds.keys():
                        number_of_rbh += 1
                weight = number_of_rbh/smallest_record_protein_count
                if weight > similarity_filter:
                    edge_map[record_accession][other_record_accession] = {'weight' : weight}
        
    G = nx.Graph(edge_map)
    
    # #clique approach
    # cliques = list(nx.find_cliques(G))
    # clique_counts = {}
    # for node in G.nodes:
    #     clique_counts[node] = len([i for i in cliques if node in i])
    # nodes_to_keep = []
    # nodes_accounted_for = set()
    # for clique in cliques:
        
    #     counts = sorted([(n, clique_counts[n]) for n in clique],
    #                     key = lambda x: x[1])
        
    #     min_count = counts[0][1]
    #     if min_count > 1:
    #         novel_nodes = [i for i in clique if i not in nodes_accounted_for]
    #         if novel_nodes == []:
    #             continue
    #         else:
    #             counts = sorted([(n, clique_counts[n]) for n in novel_nodes],
    #                     key = lambda x: x[1])
    #     nodes_accounted_for.update([i[0] for i in counts])
    #     for node, count in counts:
    #         if node not in nodes_to_keep:
    #             nodes_to_keep += [node]
    #             break
    #prune_approach
    nodes = [i for i in G.nodes if i != anchor_node]
    max_degree = 1 #does not matter, just above 0
    while True:
        temp_G = G.subgraph(nodes)
        temp_nodes = []
        degrees = []
        for node, degree in temp_G.degree:
            temp_nodes += [node]
            degrees += [degree]
        max_degree = max(degrees)
        if max_degree == 0:
            break
        delete_index = degrees.index(max_degree)#ONLY TAKE ONE NODE, NOT ALL NODES THAT HAVE SAME DEGREE
        nodes = [n for i, n in enumerate(temp_nodes) if i != delete_index]
    nodes_to_keep_pruned = nodes + [anchor_node]
    for index, acc in enumerate(nodes_to_keep_pruned):
        #write to file for clinker #
        record = cut_records[acc]
        path = f'{folder}/{record.annotations["organism"]}_{index}'
        if os.path.exists(path + '.gbk'):
            path += str(index)
        path += '.gbk'
        with open(path, 'w') as handle:
            SeqIO.write (record, 
                          handle, 
                          "genbank")
            assert os.path.isfile(path)
    plot_graph(G, f'{folder}/whole_graph.html')
    #G_clique=G.subgraph(nodes_to_keep)
    #plot_graph(G_clique, f'{folder}/clique_graph.html')
    G_prune=G.subgraph(nodes_to_keep_pruned)
    for n, d in G_prune.degree:
        if d>0:
            print (n, d)
    plot_graph(G_prune, f'{folder}/prune_graph.html')
#raise ValueError
    
'''    
clinker_accessions = []


drop = set()

for accession, other_accessions in similarity_map.items():
    #if accession in drop:
    #    print ('in drop check')
    if accession not in drop:
        for other_accession, shared_percent in other_accessions.items():
            #if other_accession in drop:
            #    print ('in drop comp')
            if other_accession not in drop and other_accession != accession:
                
                if shared_percent > 0.7:
                    drop.update([other_accession])
clinker_accessions = [i for i in similarity_map.keys() if i not in drop]
    #go to the next one that is acceptable and keep it - then start loop again but compare against both proteins and so on

#go through them all
#get rbhs between each query and hit scaffold 
        #get all of the hits grouped by host accession
        accessions = set([a.hit_def.split('__')[0] for a in record.alignments])
        grouped_alignments = {} 
        for accession in accessions:
            #grouped_alignments[accession] = []
            for alignment in record.alignments:
                keep = False
                for hsp in alignment.hsps:
                    if hsp.identities/hsp.align_length >= identity_lim:
                        keep = True
                if not keep:
                    continue
                if alignment.hit_def.split('__')[0] == accession:
                    if accession not in grouped_alignments.keys():
                        grouped_alignments[accession] = []
                    grouped_alignments[accession] += [alignment]
        raw_results[record.query] = grouped_alignments
'''
# hits = []
# for q in regions:
#     for h in hits:
#         find best hits between h and q
#         if number of h homologs/nmumber of h proteins > threshold:
#             continue
#         else:
#             hits += hits
        

# print ('completed blastp')   
#     #print (len(cut_record) == stop- start)
#     #SeqIO.write(cut_record, f"C:/Users/u03132tk/.spyder-py3/boundary_dataset/{accession}_{start}_{stop}.gbk", "genbank")
#blast - find best hit, filter by identity.  do rbh - get all the hits between smash cluster and the query cluster, then do rbjh.  show 50% 70% 90% identity requirement.  no requirmen for coverage in case the annots are broken