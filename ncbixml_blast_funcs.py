# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 10:52:06 2022

@author: u03132tk
"""
from Bio.Blast import NCBIXML
import subprocess
import sys
import time
import os
#from ModuleMapper_classes import Blast_results, RBH_result, Hit
from dataclasses import dataclass
from Bio.Blast.Record import HSP
#use this in the blast functions.  better than named tuples etc as you get named attributes
@dataclass 
class Hit:
    hit_id: str
    hsp:HSP
    query_coverage:float

@dataclass
class Blast_results:    
    #no easy way to get entire hit length so cant get hit coverage https://www.biostars.org/p/321102/
    hits:list
    query_id : str
    
    

#use this in the blast functions   z
@dataclass
class RBH_result:
    forward_result : Blast_results
    reverse_result : Blast_results
    smCOG_membership : list
    homology : list


def run_blast_process(cmd:list):
    '''this runs the BLAST+ BLASTP exe using Python via the subprocess module.  
    It uses error handlng to ease diagnostics - the return codes from BLAST+ 
    do not always match those of windows so (for example) 
    successful database creation can throw an error.'''
    try:
        return subprocess.run(cmd, check=True, capture_output=True, env={'BLASTDB_LMDB_MAP_SIZE':'1000000000'})
    except FileNotFoundError as e1:
        print ('Command:  ', cmd)
        print ('FileNotFoundError:\n', e1)
        return e1
    except subprocess.CalledProcessError as e2:
        print ('Command:  ', cmd)
        print ('CalledProcessError:')
        print ('\nError object:\n',e2)
        print ('\nStd_err:\n',e2.stderr)
        print ('\nStd_out:\n',e2.stdout)
        return e2

def makeblastdb_subprocess(makeblastdb_exe_path, input_fasta_path, db_out_path):#db_type_str param removed - redundant not used
    '''this runs the BLAST+ makeblastdb exe using Python via the subprocess module.  Input fasta is converted to BLAST db.  DB outpath is where you want the databse to be placed.
    NB - the DB outpath should have the full path up until where you want the files to be made, including the name of the files (path\to\dir\name).  
    Normally, name would have a suffix (e.g. .txt), but in this case do NOT include a file suffix.  
    Several files will be made with the same name and different suffixes, and adding a suffix does not appear to ipact databsae functionality, 
    but may lead the user to expect file functionality that does not exist - 
    ie adding a .txt suffix will not make the file behave like a .txt, as its final suffix wil be different.'''
    return run_blast_process([makeblastdb_exe_path, 
                              '-in', input_fasta_path, 
                              '-out', db_out_path, 
                              '-dbtype', 'prot'])
                              #'-parse_seqids'])#,,
                              #'-hash_index'])


def blastP_subprocess(evalue_user, query_user, blastp_exe_path, results_out_path, db_user, thread_num):
    '''this runs the biopython blastp wrapper, key details:  
        -xml output format 
        -parse deflines - this means you keep fasta heads, which is essential for coordinating rbh 
        -link to wrapper description - https://biopython.org/docs/1.75/api/Bio.Blast.Applications.html#Bio.Blast.Applications.NcbiblastpCommandline.  
    NB - the DB outpath should have the full path up until where you want the files to be made, including the name of the files (path\to\dir\name).  
    Normally, name would have a suffix (e.g. .txt), but in this case do NOT include a file suffix. 
    You are referenceing a collection of db files (.pdb, .phr, .pin, .pot, .psq, .ptf, .pto) not a single file.'''
    
    return run_blast_process([blastp_exe_path, 
                              '-out', fr'{results_out_path}', 
                              '-query', fr'{query_user}', 
                              '-db', fr'{db_user}', 
                              '-evalue', fr'{evalue_user}', 
                              '-outfmt', '5', 
                              '-num_threads', fr'{thread_num}',
                              #'-parse_deflines',#?
                              #https://microbiome.wordpress.com/research/orthologs/
                              #'-seg', 'yes',
                              #'-soft_masking', 'true',
                              #'-use_sw_tback',
                              #biopython defaults 
                              '-num_alignments', '200'])



def check_tuple_overlap(t1, t2):
    '''

    Parameters
    ----------
    t1 : tuple 
        Start, End - position indexes e.g. of a string.
    t2 : tuple
        Start, End - position indexes e.g. of a string.

    Returns
    -------
    Boolean
        Overlap between t1 and t2 (i.e. is the start and/or end one tuple within the span of the other tuple). Ignores identical start or end if other side is not overlapping ie adjacent hsps.

    '''
    t1_start_in_bool = t2[0]<t1[0]<t2[1]
    t1_end_in_bool = t2[0]<t1[1]<t2[1]
    t1_engulf_t2 = t2[0]>=t1[0] and t2[1]<=t1[1]
    return t1_start_in_bool or t1_end_in_bool or t1_engulf_t2

def find_q_coverage(hsp, query_length):
    #TODO should be alignment length
    return 100*((hsp.query_end - hsp.query_start)/query_length)



def filter_hsps(hsps : list, query_length,  min_blastscore, min_coverage):
    filtered_hsps = []
    for hsp in hsps:
        if hsp.score>= min_blastscore:
            coverage = find_q_coverage(hsp, query_length)
            if coverage >= min_coverage:
                filtered_hsps += [hsp]
    return filtered_hsps        



def collect_alignments(record, min_blastscore, min_coverage, id_field = 'hit_def'):
    hit_hsps = {}
    for alignment in record.alignments:
        if id_field == 'hit_def':
            hit_id = alignment.hit_def#.accession is missing decimels #hit_id can have gb| prefix, not sure why - and only seems to 
        elif id_field == 'hit_id':
            hit_id = alignment.hit_id
        else:
            sys.exit(f'chosen alignment id field is not recognised - {id_field}')
        filtered_hsps = filter_hsps(alignment.hsps, 
                                    record.query_length, 
                                    min_blastscore, 
                                    min_coverage)
        if len(filtered_hsps) > 0:
            if hit_id not in hit_hsps.keys():    
                hit_hsps[hit_id] = []
            hit_hsps[hit_id] += filtered_hsps
    return hit_hsps


def worst_hsp(hsp_check, hsp_compare, query_length):
    worst_hsp = 'hsp_check'
    if hsp_check.score > hsp_compare.score:
        worst_hsp = 'hsp_compare'
    elif hsp_check.score == hsp_compare.score:
        if hsp_check.expect < hsp_compare.expect:
            worst_hsp = 'hsp_compare'
        elif hsp_check.expect == hsp_compare.expect:
            check_coverage = find_q_coverage(hsp_check, query_length)
            compare_coverage = find_q_coverage(hsp_compare, query_length)
            if check_coverage > compare_coverage:
                worst_hsp = 'hsp_compare'
            #else they are equivalent, just say check is the worst
    return worst_hsp

def no_overlapping_hsps(hsps, query_length, check_query_overlap, check_hit_overlap):
    #print ('\n\n')
    assert check_hit_overlap or check_query_overlap
    removed_hsps = []
    for check_index, hsp_check in enumerate(hsps):
        #print (f'\n---CHECK {check_index}---')
        
        for compare_index, hsp_compare in enumerate(hsps):
            
            
            if check_index != compare_index:
                #if compare_index in removed_hsps:
                    #print (f'compare {compare_index} in removed_hsps')
                    
                #if check_index in removed_hsps:
                    #print (f'check {check_index} in removed_hsps')
                if   compare_index or  check_index in removed_hsps:
                    continue
                
                
                overlap = False
                if check_query_overlap:
                    overlap =  check_tuple_overlap((hsp_check.query_start, hsp_check.query_end), 
                                                   (hsp_compare.query_start, hsp_compare.query_end))
                    
                if check_hit_overlap and overlap == False:
                    #if its not false you already found it on the query
                    overlap =  check_tuple_overlap((hsp_check.sbjct_start, hsp_check.sbjct_end), 
                                                   (hsp_compare.sbjct_start, hsp_compare.sbjct_end))
                #if not overlap:
                    #print (f'{check_index} does not overlap with {compare_index}')
                if overlap:
                    #compare each hsp to every other hsp and if it overlaps, pic the best one.
                    #you then need to remove that hsp from the pool of hsps being compared, 
                    #because even if it overlaps with a different hsp it is not being considered anymore.
                    remove_hsp =  worst_hsp(hsp_check, hsp_compare, query_length)
                    if remove_hsp == 'hsp_check':
                        removed_hsps += [check_index]
                        #print (f'adding check {check_index}')
                    elif remove_hsp == 'hcp_compare':
                        removed_hsps += [compare_index]
                        #print (f'adding compare {compare_index}')
        #print (removed_hsps)
    
    assert len(set(removed_hsps)) == len(removed_hsps), f'{removed_hsps}\n{removed_hsps}\n{list(range(len(hsps)))}'                            
    return [hsp for index, hsp in enumerate(hsps) if index not in removed_hsps]  
  
def best_hsp(hsps, query_length, hit_id):
    assert len(hsps) > 0, 'you are trying to parse an empty hsp list'
    blastscores = [h.score for h in hsps]
    max_blastscore = max(blastscores)
    top_score_hsps = [hsp for hsp in hsps if hsp.score == max_blastscore]
    
    best_eval = min([hsp.expect for hsp in top_score_hsps])
    top_eval_hsps = [hsp for hsp in top_score_hsps if hsp.expect == best_eval]
    
    coverages = [find_q_coverage(h, query_length) for h in top_eval_hsps]#wont catch duplicate coverages if that occurs - unlikely but add logic for this
    query_coverage = max(coverages)
    #you dont need more than one hsp as you only want the best
    best_hsp = top_eval_hsps[coverages.index(query_coverage)]
    
    return Hit(hit_id, 
               best_hsp,
               query_coverage)
    
    

def process_hit(hit_id, hit_hsps, query_length, check_query_overlap, check_hit_overlap, best_bit_score_only):
    if check_hit_overlap or check_query_overlap:
        check_hsps = no_overlapping_hsps(hit_hsps, 
                                         query_length, 
                                         check_query_overlap, 
                                         check_hit_overlap)       
    else:
        check_hsps = hit_hsps    
    return best_hsp(check_hsps, query_length, hit_id)              
    # #print (len(hit_hsps), len(check_hsps))
    # if best_bit_score_only:
    #     hit = best_hsp(check_hsps, query_length, hit_id)
    # else:
    #     #this is potentially causing error
    #     hit = (hit_id, 
    #            sum([hsp.score for hsp in check_hsps]), 
    #            min([hsp.expect for hsp in check_hsps]),#summing evals isnt the riht way to go about this - maybe mean is better
    #            sum([find_q_coverage(hsp, query_length) for hsp in check_hsps])
    #            )
    # return hit

def best_hit(hits, query_id):
    if len(hits) == 0:
        max_coverage_hits = []
    else:
        #hits = [(hit_id, hit_hsp, query_length), ()...]
        hit_scores = [hit.hsp.score for hit in hits]
        max_score = max(hit_scores)
        max_score_hits = [hit for hit in hits if hit.hsp.score == max_score]
        #print (f'\nthere are {len(max_score_hits)} hits')
        hit_evals = [hit.hsp.expect for hit in max_score_hits]
        min_eval = min(hit_evals)
        min_eval_hits = [hit for hit in max_score_hits if hit.hsp.expect == min_eval]
        #print (f'there are {len(min_eval_hits)} hits')
        hit_coverages = [hit.query_coverage for hit in min_eval_hits]
        max_coverage = max(hit_coverages)
        #print (max_coverage)
        #print ([(hit, hit[1].query_end, hit[1].query_start) for hit in min_eval_hits])
        #max_coverage_hit = min_eval_hits[hit_coverages.index(max_coverage)]#pick one, will inore duplicates here
        #you could get more than one hit
        max_coverage_hits = [hit for hit  in min_eval_hits if hit.query_coverage == max_coverage]
        '''
        make this a list and append any obj with same max coverage - needed for multicopy proteins.  
        could add a lambda expression for these params - e.g pick any within 1% of top hit
        make return a dict obj, with one hit per query
        '''
        #except TypeError as e:
            #print (type(min_eval_hits), type(hit_coverages), type(max_coverage), type(hit_coverages.index(max_coverage)))
            #raise e
    #print (f'there are {len(max_coverage_hits)} hits')
        
    return Blast_results (hits = max_coverage_hits,
                          query_id = query_id)
                               


def parse_xml(fw_xml_results_path, best_hit_only,
              check_hit_overlap, check_query_overlap, min_blastscore, min_coverage):
    queries = []
    
    with open(fw_xml_results_path, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        #for each query r, go through each hit
        for r in blast_records:
            #collect all alignments, extract hsps, map hsp lists to their respective hits
            hit_alignments = collect_alignments(r, 
                                                min_blastscore, 
                                                min_coverage)
            
            #go through the hsps for each hit and process to a final set of hit statsictics
            hits = []
            for hit_id, hit_hsps in hit_alignments.items():#r.alignments:
                #print (hit_id, 'has', len(hit_hsps), 'hsps')
                # if 'RS55105' in r.query_id or 'SMCOG10183__YP_003112868.1' in r.query_id or "M878_RS80075" in r.query_id:
                #     if 'RS55105' in hit_id or 'SMCOG10183__YP_003112868.1' in hit_id or 'M878_RS80075' in hit_id:
                #         print (f'\nQUERY {r.query_id} has {len(hit_hsps)} for hit {hit_id} before processing')
                         
                hits += [process_hit(hit_id, 
                                     hit_hsps, 
                                     r.query_length, 
                                     check_query_overlap, 
                                     check_hit_overlap, 
                                     best_hit_only)
                         ]
                # if 'RS55105' in r.query_id or 'SMCOG10183__YP_003112868.1' in r.query_id or "M878_RS80075" in r.query_id:
                #     if 'RS55105' in hit_id or 'SMCOG10183__YP_003112868.1' in hit_id or 'M878_RS80075' in hit_id:
                #         print (f'QUERY {r.query_id} has {len(hit_hsps)} for hit {hit_id} after processing')
                #         for hsp in hit_hsps:
                #             print (f'score {hsp.score} eval {hsp.expect} from query {hsp.query_start} to {hsp.query_end} and hit {hsp.sbjct_start} to {hsp.sbjct_end}')
                #     # elif 'SMCOG10183__YP_003112868' in r.query_id:
                #     #     print (f'QUERY {r.query_id} has {len(hit_hsps)} for hit {hit_id} after processing - not normal loop ')
                #     #     for hsp in hit_hsps:
                #     #         print (f'score {hsp.score} eval {hsp.expect} from query {hsp.query_start} to {hsp.query_end} and hit {hsp.sbjct_start} to {hsp.sbjct_end}')
            #def parse hit list
            
            #TODO  - needs to handle no hit
            queries.append(best_hit(hits, 
                                    r.query))
    return queries

#TODO if this is too slow move list comprehensions to loops.  will save repeated looping but will cost readability.


def makedb_blastp_and_parse(makeblastdb_exe_path,
                            blastp_exe_path,
                            input_db_fasta_path,
                            db_outpath,
                            query_user,
                            results_out_path,
                            thread_num,
                            evalue_user,
                            best_hit_only,
                            min_query_perc_coverage, 
                            min_blastscore,
                            check_query_overlap,
                            check_hit_overlap):
    start_db=time.time()
    makeblastdb_subprocess(makeblastdb_exe_path, 
                           input_db_fasta_path, 
                           db_outpath)
    print (f'make db:  {time.time()-start_db}')
    start_blastp=time.time()
    blastP_subprocess (evalue_user, 
                       query_user, 
                       blastp_exe_path, 
                       results_out_path, 
                       db_outpath, 
                       thread_num)
    print (f'complete blastp:  {time.time()-start_blastp}')
    start_parse=time.time()
    best_hits=parse_xml(results_out_path, 
                        best_hit_only, 
                        check_hit_overlap, 
                        check_query_overlap, 
                        min_blastscore, 
                        min_query_perc_coverage)
    print (f'complete parse:  {time.time()-start_parse}')
    return best_hits

#ADD propertyword_size
def blastp_and_parse(blastp_exe_path, 
                     db_path, 
                     query_user, 
                     results_out_path,
                     thread_num,
                     evalue_user,
                     best_hit_only, 
                     check_hit_overlap, 
                     check_query_overlap, 
                     min_blastscore, 
                     min_query_perc_coverage):
    '''
    

    Parameters
    ----------
    error_file : TYPE
        DESCRIPTION.
    blastp_exe_path : TYPE
        DESCRIPTION.
    db_path : TYPE
        DESCRIPTION.
    query_user : TYPE
        DESCRIPTION.
    results_out_path : TYPE
        DESCRIPTION.
    thread_num : TYPE, optional
        DESCRIPTION. The default is os.cpu_count().
    evalue_user : TYPE, optional
        DESCRIPTION. The default is 1*(10**-6).
    best_hit_only : TYPE, optional
        DESCRIPTION. The default is False.
    min_hit_perc_coverage : TYPE, optional
        DESCRIPTION. The default is 0.
    min_query_perc_coverage : TYPE, optional
        DESCRIPTION. The default is 0.
    check_overlap : TYPE, optional
        DESCRIPTION. The default is True.
    overlap_hit_check_only : TYPE, optional
        DESCRIPTION. The default is False.
    overlap_query_check_only : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    best_hits : TYPE
        DESCRIPTION.

    '''
    start_blastp=time.time()
    blastP_subprocess (evalue_user, query_user, blastp_exe_path, results_out_path, db_path, thread_num)
    print (f'complete blastp:  {time.time()-start_blastp}')
    start_parse=time.time()
    best_hits=parse_xml(results_out_path, 
                        best_hit_only, 
                        check_hit_overlap, 
                        check_query_overlap, 
                        min_blastscore, 
                        min_query_perc_coverage)
    print (f'complete parse:  {time.time()-start_parse}')
    return best_hits

def two_way_compare_results(results1,results2):
    different_query_hits=[]
    for query_t,query_b in zip(results1, results2):
        if query_t!=query_b:
            different_query_hits.append([query_t, query_b])
    smcog_mismatch=[]
    for i in different_query_hits:
        if i[0][0]!=i[1][0]:
            smcog_mismatch.append(i)
    return different_query_hits,smcog_mismatch


def check_fasta(fasta_old, fasta_new):
    '''
    checks fasta deflines is associated with same seq in both fastas - does not check order
    '''
    with open(fasta_old,'r') as file_in:
        original_fasta_data=[line.strip() for line in file_in.readlines()]
        original_deflines = [i for i in original_fasta_data if '>' in i]
        duplicate_deflines_o = len(set(original_deflines))> len(original_deflines)
        if duplicate_deflines_o:
            sys.exit('duplicate deflines in forward database')
        with open(fasta_new,'r') as file_out:
            new_fasta_data=[line.strip() for line in file_out.readlines()]
            new_deflines = [i for i in new_fasta_data if '>' in i]
            duplicate_deflines_n=len (set(new_deflines))> len(new_deflines)
            if duplicate_deflines_n:
                sys.exit('duplicate deflines in reverse query')
            for defline in new_deflines:
                old_index=original_fasta_data.index(defline)+1
                new_index=new_fasta_data.index(defline)+1
                sequence_new=new_fasta_data[new_index]
                sequence_old=original_fasta_data[old_index]
                if sequence_new!=sequence_old:
                    return [False, (defline, sequence_new, sequence_old)]
    return [True]
    

def build_reciprical_blastp_query(original_db_fasta_path, new_query_fasta_path, parsed_xml_list):
    start_write=time.time()
    already_written=[]
    with open(original_db_fasta_path,'r') as file_in:
        original_fasta_data=[line.strip() for line in file_in.readlines()]
    with open(new_query_fasta_path,'w') as file_out:
        for result in parsed_xml_list:
            if len(result.hits)>0:
                for hit in result.hits:
                    defline = f'>{hit.hit_id}'
                    if defline not in already_written:
                        sequence_index=original_fasta_data.index(defline)+1
                        sequence=original_fasta_data[sequence_index]
                        file_out.write(f'{defline}\n{sequence}\n')
                        already_written.append(defline)
    print (f'write reciprical query:  {time.time()-start_write}')
    start_check=time.time()
    correct_fasta_transfer=check_fasta(original_db_fasta_path, new_query_fasta_path)
    print (f'check new query of database hits against original database fasta:  {time.time()-start_check}')
    print (f'new query correct:  {correct_fasta_transfer}')
    return correct_fasta_transfer

     
#potentially look here for additional blast filtering settings - not they use 10-5 not 10-6 evalue!
#https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3196566/ 

# you are loosing smcog assignment for E4U91_RS34740 in SZN bgc 29 extended 5kb either side - make sure you extend so you dont overwrite files.
# on further analysis it looks like the reverse hit for the smcog is not the same as the query.
# find_q_coverage is giveing 409% coverage, implying q start - end is 4x larger than query
# i think error lies in hit processing - make sure your hsps are definitely not overlapping.  

# it is!  best hit only replicates original results  


def blast_rbh_optional_reverse(fw_query_fasta,
                               fw_db_fasta_path,
                               fw_db_path,
                               fw_results_out_path,
                               rev_query_fasta_path,
                               rev_db_fasta_path,
                               rev_db_outpath,
                               rev_results_out_path, 
                               makeblastdb_exe_path,
                               blastp_exe_path, 
                               best_hit_only=True, 
                               min_query_perc_coverage=0, 
                               check_hit_overlap = True, 
                               check_query_overlap = True, 
                               min_blastscore = 0,
                               thread_num=os.cpu_count(),
                               evalue_user= 1*(10-5), 
                               allowed_homology = 'ortholog'):#add all the blast options 
    start_f=time.time()  
    forward_hits=blastp_and_parse(blastp_exe_path, 
                                  fw_db_path, 
                                  fw_query_fasta, 
                                  fw_results_out_path,
                                  thread_num,
                                  evalue_user,
                                  best_hit_only, 
                                  check_hit_overlap, 
                                  check_query_overlap, 
                                  min_blastscore, 
                                  min_query_perc_coverage)
    print (f'forward blastp total time\t{time.time()-start_f} seconds')
    #reverse_query_file=build_reciprical_blastp_query(fw_db_fasta_path, rev_query_fasta_path, forward_hits)
    start_r=time.time() 
    reverse_query_file=build_reciprical_blastp_query(fw_db_fasta_path, rev_query_fasta_path, forward_hits)
    if reverse_query_file[0]==True:
        #reverse_db_fasta_path=query_user
        reverse_hits=makedb_blastp_and_parse(makeblastdb_exe_path,
                                             blastp_exe_path,
                                             rev_db_fasta_path,
                                             rev_db_outpath,
                                             rev_query_fasta_path,
                                             rev_results_out_path,
                                             thread_num,
                                             evalue_user,
                                             best_hit_only,
                                             min_query_perc_coverage, 
                                             min_blastscore,
                                             check_query_overlap,
                                             check_hit_overlap)
        print (f'reverse blastp total time\t{time.time()-start_r} seconds')
    else:
        sys.exit('error in reverse query file') 
    #error here in second json
    start_table=time.time()
    evo_table={}
    
    #seems rare but you can have a fw hit that does not have any reverse hits when used as a query - e.g. SMCOG19798__YP_001159858.1 hit of H7H31_RS30715 for NZ_CP060404.1_r42.json
                # if bgc_db_defline_f == None:
                #     continue
    
    for f_result in forward_hits:
        evo_table[f_result.query_id] = []
        if len(f_result.hits)>0:#none = no hits
            at_least_one_ortholog = False
            if 'ortholog' in allowed_homology:
                for f_hit in f_result.hits:
                    for r_result in reverse_hits:
                        if r_result.query_id == f_hit.hit_id:
                            for r_hit in r_result.hits:
                                if f_result.query_id == r_hit.hit_id:
                                    evo_table[f_result.query_id].append(RBH_result(#f_query = f_result.query_id, 
                                                                                   forward_result = f_hit, 
                                                                                   reverse_result = r_hit, 
                                                                                   smCOG_membership = [],
                                                                                   homology = 'ortholog'))
                                    at_least_one_ortholog = True 
            #if there is an ortholog assigned do not look for homologs even if they are allowed
            #this means you need to check all of the hits for othologs first 
            if 'homolog' in allowed_homology:        
                if not at_least_one_ortholog:
                    #these are equally acceptable hits
                    for f_hit in f_result.hits:  
                        for r_result in reverse_hits:
                            #find the hit of the reverse query 
                            if f_hit.hit_id == r_result.query_id:
                                for r_hit in r_result.hits:
                                    evo_table[f_result.query_id].append(RBH_result(#f_query = f_result.query_id, 
                                                                                   forward_result = f_hit, 
                                                                                   reverse_result = r_hit, 
                                                                                   smCOG_membership = [],
                                                                                   homology = 'homolog'))
    print (f'table makeup took {time.time()-start_table} seconds')
    return evo_table

'''
evo_table=[]
    for hit_f in forward_hits:
        bgc_query_defline_f = hit_f[0]
        smcog_db_defline_f=hit_f[1]

        if smcog_db_defline_f != None:#none = no hits
            split_defline_f_db=smcog_db_defline_f.split('__')
            smcog_id=split_defline_f_db[0]
            relationship_assigned=False
            #if 'RS55105' in bgc_query_defline_f or "M878_RS80075" in bgc_query_defline_f:
            #    print (hit_f)
                    
            for hit_r in reverse_hits:
                #if 'SMCOG10183__YP_003112868.1' in smcog_db_defline_f:
                #    print (hit_r)
                smcog_query_defline_r=hit_r[0]
                bgc_db_defline_f=hit_r[1]
                #this accounts for same protein as well, not just same smcog
                check_f = bgc_query_defline_f.replace('dbj|', '').replace('gb|', '').replace('|', '') == bgc_db_defline_f.replace('dbj|', '').replace('gb|', '').replace('|', '')
                check_r = smcog_db_defline_f.replace('dbj|', '').replace('gb|', '').replace('|', '') == smcog_query_defline_r.replace('dbj|', '').replace('gb|', '').replace('|', '')
                if check_f and check_r:
                    #evo_table.append([bgc_id,bgc_orf,smcog_id,'ortholog'])
                    evo_table.append([bgc_query_defline_f,smcog_id,'ortholog'])
                    relationship_assigned=True
                    break
                    #this condition can only be satisfied once as a query can only have one hit and there should be no query duplicates.  
                    #although there can be duplicate smcog hits (i.e. multiple smcog_hit==smcog_query), 
                    #there will only be one instance where the dem query is matched (dem_query==dem_hit).  
                    #This means that if there are multiple hits to a single DEM protein from the smcog queries, only one of these will be the 
                    #smcof that the original query matched in the forward blastp. 
            # you get all the info you need here
            # you will need to find the query coverage here, as well as alignment 
            if not relationship_assigned:
                evo_table.append([bgc_query_defline_f,smcog_id,'homolog'])
    print (f'table makeup took {time.time()-start_table} seconds')
    return evo_table
'''


#if you cant get this working straight away go back to the old db that was commented out - but make it a priority to get this owrking

if __name__ == '__main__':
    f_xml_path = r'C:/Users/u03132tk/.spyder-py3/ModuleMapper/Backend/Intermediate_files/forward_blast/BLASTP_results/NZ_BIFH00000000/json_NZ_BIFH00000000__genome_number_1__genome_id_NZ_BIFH01000013__for_fw_results.xml'
    with open(f_xml_path, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)
    f = parse_xml(f_xml_path, True, True,True,0,0 )
    with open(f_xml_path, 'r') as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
    #sys.exit()
    
    
    r_xml_path = r'C:/Users/u03132tk/.spyder-py3/ModuleMapper/Backend/Intermediate_files/reverse_blast/BLASTP_results/NZ_BIFH00000000/json_NZ_BIFH00000000__genome_number_1__genome_id_NZ_BIFH01000013__for__results_reverse.xml'
    r = parse_xml(r_xml_path, True, True,True,0,0 )
    evo_table={}
    for f_result in f:
        evo_table[f_result.query_id] = []
        if len(f_result.hits)>0:#none = no hits
            at_least_one_ortholog = False
            for f_hit in f_result.hits:
                f_hit_id = f_hit.hit_id
                for r_result in r:
                    if r_result.query_id == f_hit.hit_id:
                        for r_hit in r_result.hits:
                            if f_result.query_id == r_hit.hit_id:
                                evo_table[f_result.query_id].append(RBH_result(#f_query = f_result.query_id, 
                                                                               forward_result = f_hit, 
                                                                               reverse_result = r_hit, 
                                                                               smCOG_membership = [],
                                                                               homology = 'ortholog'))
                                at_least_one_ortholog = True 
                    
            if not at_least_one_ortholog:
                #these are equally acceptable hits
                for f_hit in f_result.hits:  
                    f_hit_id = f_hit.hit_id
                    f_hits = f_hit
                    for r_result in r:
                        r_query = r_result.query_id
                        #find the hit of the reverse query 
                        if f_hit_id == r_query:
                            for r_hit in r_result.hits:
                                evo_table[f_result.query_id].append(RBH_result(#f_query = f_result.query_id, 
                                                                               forward_result = f_hit, 
                                                                               reverse_result = r_hit, 
                                                                               smCOG_membership = [],
                                                                               homology = 'homolog'))
               