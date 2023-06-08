import argparse
from qnorm_fast import qnorm_reducesort
from qnorm_fast import qnorm
import scipy.sparse
import numpy as np
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('merge_0min', type=str, help="merge 0 min iced counts")
    parser.add_argument('merge_25min', type=str, help="merge 25 min iced counts")
    parser.add_argument('merge_60min', type=str, help="merge 60 min iced counts")
    parser.add_argument('merge_120min', type=str, help="merge 120 min iced counts")
    parser.add_argument('merge_240min', type=str, help="merge 240 min iced counts")

    parser.add_argument('chr', type=str, help = "chromosome")
    parser.add_argument('resolution',type = int, help = 'resolution')
    parser.add_argument('threshold',type = float, help = 'threshold')

    args = parser.parse_args()

    command_removal = "rm " +  args.merge_0min[:-4] + "_qnormed*_merge*" + args.chr + "*.npz"
    os.system(command_removal)
    command_removal = "rm " +  args.merge_25min[:-4] + "_qnormed*_merge*" + args.chr + "*.npz"
    os.system(command_removal)
    command_removal = "rm " +  args.merge_60min[:-4] + "_qnormed*_merge*" + args.chr + "*.npz"
    os.system(command_removal)
    command_removal = "rm " +  args.merge_120min[:-4] + "_qnormed*_merge*" + args.chr + "*.npz"
    os.system(command_removal)
    command_removal = "rm " +  args.merge_240min[:-4] + "_qnormed*_merge*" + args.chr + "*.npz"
    os.system(command_removal)


    merge_0min = scipy.sparse.load_npz(args.merge_0min)
    merge_25min = scipy.sparse.load_npz(args.merge_25min)
    merge_60min = scipy.sparse.load_npz(args.merge_60min)
    merge_120min = scipy.sparse.load_npz(args.merge_120min)
    merge_240min = scipy.sparse.load_npz(args.merge_240min) 

    maximum_size = np.max([merge_0min.shape[0],merge_25min.shape[0],merge_60min.shape[0],merge_120min.shape[0],merge_240min.shape[0]])

    merge_0min._shape = (maximum_size,maximum_size)
    merge_25min._shape = (maximum_size,maximum_size)
    merge_60min._shape = (maximum_size,maximum_size)
    merge_120min._shape = (maximum_size,maximum_size)
    merge_240min._shape = (maximum_size,maximum_size)

    merge_0min.tocoo()
    merge_25min.tocoo()
    merge_60min.tocoo()
    merge_120min.tocoo()
    merge_240min.tocoo()

    merge_0min_coo_coord = zip(merge_0min.row,merge_0min.col)
    merge_25min_coo_coord = zip(merge_25min.row,merge_25min.col)
    merge_60min_coo_coord = zip(merge_60min.row,merge_60min.col)
    merge_120min_coo_coord = zip(merge_120min.row,merge_120min.col)
    merge_240min_coo_coord = zip(merge_240min.row,merge_240min.col)


    #collect total (row,col) coordinates with values found across all replicates
    coordinates = list(set(merge_0min_coo_coord).union(merge_25min_coo_coord,merge_60min_coo_coord,merge_120min_coo_coord,merge_240min_coo_coord))
    coordinates = np.array(coordinates) 
    
    ext_data_rep = np.zeros((len(coordinates),5))

    merge_0min_csr = merge_0min.tocsr()
    merge_25min_csr = merge_25min.tocsr()
    merge_60min_csr = merge_60min.tocsr()
    merge_120min_csr = merge_120min.tocsr()
    merge_240min_csr = merge_240min.tocsr()

    
    #store values per replicate at corresponding sorted coordinate (rows = total_pixels, cols = replicates)   
    ext_data_rep[:,0] = merge_0min_csr[coordinates[:,0],coordinates[:,1]]
    ext_data_rep[:,1] = merge_25min_csr[coordinates[:,0],coordinates[:,1]]
    ext_data_rep[:,2] = merge_60min_csr[coordinates[:,0],coordinates[:,1]]
    ext_data_rep[:,3] = merge_120min_csr[coordinates[:,0],coordinates[:,1]]
    ext_data_rep[:,4] = merge_240min_csr[coordinates[:,0],coordinates[:,1]]

    print "total coordinates: ", len(coordinates)

    np.set_printoptions(suppress=True) 

    #convert nans to zero for nan safety (avoid affecting rank sorting for qnorming)
    ext_data_rep[np.isnan(ext_data_rep)] = 0 

    #remove all rows <= threshold percent of nonzeros in row 
    row_removal = np.where(np.count_nonzero(ext_data_rep, axis=1) <= args.threshold*ext_data_rep.shape[1])

    if len(row_removal) > 0:
        row_removal = row_removal[0]
        print "rows to be removed: ", row_removal
        print "number of rows: ", len(row_removal)
        ext_data_rep = np.delete(ext_data_rep,row_removal,axis=0)
        coordinates = np.delete(coordinates,row_removal,axis=0).tolist()
    else:
        coordinates = coordinates.tolist()


    ref_index = None
    tie='min'  #lowest rank
 
    # qnorm matrices
    qnormed_matrices = qnorm_reducesort(ext_data_rep,tie=tie)
    qnormed_matrices[~np.isfinite(qnormed_matrices)] = np.nan


    rows = np.array(list(zip(*coordinates)[0]))
    cols = np.array(list(zip(*coordinates)[1]))

    print "rows length = ", len(rows)
    print "cols length = ", len(cols)
    print "qnormed_matrices length = ", len(qnormed_matrices[0])

    qnorm_merge_0min_triu = scipy.sparse.coo_matrix((qnormed_matrices[:,0],(rows,cols)),shape=merge_0min._shape).tocsr()
    qnorm_merge_25min_triu = scipy.sparse.coo_matrix((qnormed_matrices[:,1],(rows,cols)),shape=merge_25min._shape).tocsr()
    qnorm_merge_60min_triu = scipy.sparse.coo_matrix((qnormed_matrices[:,2],(rows,cols)),shape=merge_60min._shape).tocsr()
    qnorm_merge_120min_triu = scipy.sparse.coo_matrix((qnormed_matrices[:,3],(rows,cols)),shape=merge_120min._shape).tocsr()
    qnorm_merge_240min_triu = scipy.sparse.coo_matrix((qnormed_matrices[:,4],(rows,cols)),shape=merge_240min._shape).tocsr()    


    #print out  matrix
    scipy.sparse.save_npz(args.merge_0min[:-4] + "_qnormed_merge_allcond_" + str(args.threshold) + '_' + str(args.resolution) + '_0min_'  + args.chr + "_below0.5.npz",qnorm_merge_0min_triu)
    scipy.sparse.save_npz(args.merge_25min[:-4] + "_qnormed_merge_allcond_" + str(args.threshold) + '_' + str(args.resolution) + '_25min_' + args.chr + "_below0.5.npz",qnorm_merge_25min_triu)
    scipy.sparse.save_npz(args.merge_60min[:-4] + "_qnormed_merge_allcond_" + str(args.threshold) + '_' +  str(args.resolution) + '_60min_' + args.chr + "_below0.5.npz",qnorm_merge_60min_triu)
    scipy.sparse.save_npz(args.merge_120min[:-4] + "_qnormed_merge_allcond_" + str(args.threshold) + '_' +  str(args.resolution) + '_120min_' + args.chr + "_below0.5.npz",qnorm_merge_120min_triu)
    scipy.sparse.save_npz(args.merge_240min[:-4] + "_qnormed_merge_allcond_" + str(args.threshold) + '_' +  str(args.resolution) + '_240min_' + args.chr + "_below0.5.npz",qnorm_merge_240min_triu)

if __name__ == "__main__":
    main()
