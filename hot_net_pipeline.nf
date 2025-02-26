

params.network_dir="$projectDir/data/networks/"
params.hot_net_dir="$projectDir/"
params.out = "$projectDir/out/"
params.ms = 1000
params.nb_permut = 100
params.lsb = 10


process Build_sim_matrix{
/*
    This process build the similarity matrix

    param network_edges: path to a file describing the network edge
    out sim_matrix: path tot he generated similarity matrix in h5 format
    out path tot he beta file
*/
    publishDir params.out+"/intermediate", mode: 'copy'
    
    input:
        path network_edges

    output:
        path "*.sim_matrix.h5" , emit: sim_matrix
        path "*.beta.txt", emit: beta


    script:
    """
    python $params.hot_net_dir/src/construct_similarity_matrix.py -i $network_edges  -o ${network_edges}.sim_matrix.h5 -bof ${network_edges}.beta.txt
    """
}

process Find_permutation_bin{
/*
    This process find the permutaiton bin for th epermutaiton analysis

    param score: path to the score file
    param index_gene: path to the gene index file
    param edge_list: path to the edge list file
    out: score bin file
*/
    publishDir params.out+"/intermediate", mode: 'copy'

    input:
        path score
        path index_gene
        path edge_list

    output:
        path "*.score_bins.tsv", emit: score_bins

    script:
    """
    python $params.hot_net_dir/src/find_permutation_bins.py  -gsf $score -igf $index_gene -elf $edge_list -ms $params.ms -o ${edge_list}_${score}.score_bins.tsv
    """

}

process Create_n_suffixed_path{
/*
    This technical process is used to create n file path sufixed with an index i (1-->n)

    param base_file: path to the base file we want to index
    param n_iter: number of ieration (size of N)
    out list_file_iter: path to a file storing the n sufiixed files
*/
    publishDir params.out+"/intermediate", mode: 'copy'

    input:
        path base_file
        val n_iter
    output:
        path "*.list_files_iter", emit: list_file_iter

    """
    #! python

    with open("${base_file}.list_files_iter", "w") as file_handler:
        for i in range(1,${n_iter}+1):
            file_handler.write("${base_file}_iter_"+str(i)+".tsv\\n")
    
    """
}

process Permutation_score{
/*
    This prcoess generate a permitation of the scores

    param permut_info: list fo value [val1,val2,val3,val4] use as paramter of the script run by the process.
    out permut_score:path chanel with the indexed score
*/

    //publishDir params.out+"/intermediate/", mode: 'copy'

     input:
        val permut_info
     
     output:
        path "*_iter_*.tsv", emit: permut_score
     script:
     """
        python $params.hot_net_dir/src/permute_scores.py -i ${permut_info[1]} -bf ${permut_info[4]} -s ${permut_info[3]} -o ${permut_info[2]}
     """
}

process Build_hierachy_initial{
/*    
    This process generate the hierachy file for the non permuted data

    param sim_matrix: similarity matrix
    param index_genes: index genes 
    param scores:path to the score file
    out hier_edge_list: chanel with the hierarchy_edge_list file
    out hier_ind_gene: chanel with the hier_ind_gene file
*/
    publishDir params.out+"/intermediate/", mode: 'copy'

    input:
        path sim_matrix
        path index_genes
        path scores

    output:
        path "*-hierarchy_edge_list.tsv", emit: hier_edge_list
        path "*-hierarchy_index_gene.tsv", emit: hier_ind_gene
    """
     python $params.hot_net_dir/src/construct_hierarchy.py \
                                   -smf $sim_matrix \
                                   -igf $index_genes \
                                   -gsf $scores \
                                   -helf ${index_genes}-hierarchy_edge_list.tsv \
                                   -higf ${index_genes}-hierarchy_index_gene.tsv
    """
}

process Build_hierachies{
/*    
    This process generate the hierachy file for the  permuted data

    param input_chanel: chanel with a list fo values [val1,val2,val3,val4] used as parameter for the script
    out hier_edge_list: chanel with the hierarchy_edge_list files
    out hier_ind_gene: chanel with the hier_ind_gene files                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
*/
    
   //publishDir params.out+"/intermediate/", mode: 'copy'

    input:
        val input_chanel

    output:
        path "*-hierarchy_edge_list.tsv", emit: hier_edge_list
        path "*-hierarchy_index_gene.tsv", emit: hier_ind_gene

    script:
     """
     python $params.hot_net_dir/src/construct_hierarchy.py -smf  ${input_chanel[3]} \
                                       -igf  ${input_chanel[4]} \
                                       -gsf  ${input_chanel[1]}\
                                       -helf ${input_chanel[0]}_${input_chanel[2]}-hierarchy_edge_list.tsv \
                                       -higf ${input_chanel[0]}_${input_chanel[2]}-hierarchy_index_gene.tsv
     """

}

process Pocess_hierarchy{

    publishDir params.out+"/result/", mode: 'copy'

    input:
        path hier_edge_list
        path hier_index_gene
        path score
        path net
        val permut_edge
        val permut_i_gene

    output:
         path "clusters_*", emit: cluster
         path "sizes_*", emit: sizes

    script:
    """
     python $params.hot_net_dir/src/process_hierarchies.py -oelf $hier_edge_list -oigf $hier_index_gene -pelf $permut_edge \
                                                           -pigf $permut_i_gene  -lsb $params.lsb -cf clusters_${net}_${score}.tsv \
                                                           -pl $net $score -pf sizes_${net}_${score}.pdf 
    """
}


nextflow.enable.dsl=2

workflow{

    network_edges = Channel.fromPath(params.network_dir+'/edge_mippie.tsv')
    network_score = Channel.fromPath(params.network_dir+'/score_mippie.tsv')
    index_gene = Channel.fromPath(params.network_dir+'/mippie_gene_index.tsv')

    sim_matrix = Build_sim_matrix(network_edges)
    score_bins = Find_permutation_bin(network_score,index_gene,network_edges).collect().map{it[0]} 
    file_nb_perm = Create_n_suffixed_path(network_score,params.nb_permut)
    chanel_to_go = file_nb_perm.splitText()

    // Use combine() to pair score_bins to all elment of the chanel
    chanel_iter_perm = chanel_to_go.combine(score_bins).map { pair -> 
        def (file, bin) = pair
        [file.toString().split("_iter_")[0], file, file.toString().split("_iter_")[-1].split(".tsv")[0], bin]
    }
    network_score_id = network_score.map { [it.toString().split("/")[-1], it]}
    pair_chanel_iter_perm = network_score_id.combine(chanel_iter_perm, by: 0)
    permutations = Permutation_score(pair_chanel_iter_perm)

    sim_mat = sim_matrix.sim_matrix.collect().map{it[0]} 
    map_sim_mat_perm = permutations.combine(sim_mat).map { pair -> 
        def (file, sim_mat) = pair
        [file.toString().split("_iter_")[0].split("/")[-1], file, file.toString().split("_iter_")[-1].split(".tsv")[0], sim_mat ]
    }
    ind_gene = index_gene.collect().map{it[0]} 
    ind_gene_map = map_sim_mat_perm.combine(ind_gene).map { elem -> 
        def (val1,val2,val3,val4,ind_gene) = elem
        [val1,val2,val3,val4, ind_gene]
    }

    // build hierarchy
    hier_init = Build_hierachy_initial(sim_matrix.sim_matrix, index_gene, network_score)
 
    hier_arch_chan = Build_hierachies(ind_gene_map)

    hier_edge_lst_param = hier_arch_chan.hier_edge_list.collect().map { it.join(" ") } 
    hier_ind_gene_param = hier_arch_chan.hier_ind_gene.collect().map { it.join(" ") } 

    //edg_lst_id_ind_g_id.view()
    hier_arch_chan.hier_edge_list.view() 
    hier_arch_chan.hier_edge_list.flatten().view() 
    Pocess_hierarchy(hier_init.hier_edge_list, hier_init.hier_ind_gene,
                     network_score, network_edges, hier_edge_lst_param,
                     hier_ind_gene_param)
}




