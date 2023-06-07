import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot
import scipy.sparse
import numpy as np
import argparse
import glob
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('pathway_r1', type=str, help="subfolder rep1")
    parser.add_argument('pathway_r2', type=str, help="subfolder rep2")
    parser.add_argument('pathway_r3', type=str, help="subfolder rep3")
    parser.add_argument('pathway_r4', type=str, help="subfolder rep4")
    parser.add_argument('pathway_r5', type=str, help="subfolder rep5")
    parser.add_argument('chr', type=str, help = "chromosome")
    parser.add_argument('resolution', type=int, help = "resolution")
    parser.add_argument('type',type=str, help = "experimental type: kr,iced,kr+qnorm")
    parser.add_argument('threshold_pre',type=str, help = "cutoff for qnorm")
    parser.add_argument('threshold',type=float, help = "count cutoff for qnorming")

 
    args = parser.parse_args()

    file1 = args.pathway_r1 + '/final_matrix_rep1'
    file2 = args.pathway_r2 + '/final_matrix_rep2'
    file3 = args.pathway_r3 + '/final_matrix_rep3'
    file4 = args.pathway_r4 + '/final_matrix_rep4'
    file5 = args.pathway_r5 + '/final_matrix_rep5'

    if os.path.isfile(file1):
        os.remove(file1)
    if os.path.isfile(file2):
        os.remove(file2)
    if os.path.isfile(file3):
        os.remove(file3)
    if os.path.isfile(file4):
        os.remove(file4)
    if os.path.isfile(file5):
        os.remove(file5)

    if args.type == 'kr' or args.type == 'iced' or args.type == 'raw':
        if args.resolution == 10000: 
            find_path_r1 = args.pathway_r1 + '/' + args.chr + '_' + args.type + '.npz'
            find_path_r2 = args.pathway_r2 + '/' + args.chr + '_' + args.type + '.npz'
            find_path_r3 = args.pathway_r3 + '/' + args.chr + '_' + args.type + '.npz'
            find_path_r4 = args.pathway_r4 + '/' + args.chr + '_' + args.type + '.npz'
            find_path_r5 = args.pathway_r5 + '/' + args.chr + '_' + args.type + '.npz'
        else:
            find_path_r1 = args.pathway_r1 + '/' + args.chr + '_' + str(args.resolution) + '_' + args.type + '.npz'
            find_path_r2 = args.pathway_r2 + '/' + args.chr + '_' + str(args.resolution) + '_' + args.type + '.npz'
            find_path_r3 = args.pathway_r3 + '/' + args.chr + '_' + str(args.resolution) + '_' + args.type + '.npz'
            find_path_r4 = args.pathway_r4 + '/' + args.chr + '_' + str(args.resolution) + '_' + args.type + '.npz'
            find_path_r5 = args.pathway_r5 + '/' + args.chr + '_' + str(args.resolution) + '_' + args.type + '.npz'

        print "type is ", args.type
        print "grabbing file: ",find_path_r1
        print "grabbing file: ",find_path_r2
        print "grabbing file: ",find_path_r3
        print "grabbing file: ",find_path_r4
        print "grabbing file: ",find_path_r5

        matrix_rep1 = scipy.sparse.load_npz(find_path_r1)
        matrix_rep1 = matrix_rep1.tocsr()
        matrix_rep2 = scipy.sparse.load_npz(find_path_r2)
        matrix_rep2 = matrix_rep2.tocsr()
        matrix_rep3 = scipy.sparse.load_npz(find_path_r3)
        matrix_rep3 = matrix_rep3.tocsr()
        matrix_rep4 = scipy.sparse.load_npz(find_path_r4)
        matrix_rep4 = matrix_rep4.tocsr()
        matrix_rep5 = scipy.sparse.load_npz(find_path_r5)
        matrix_rep5 = matrix_rep5.tocsr()

    else:
        print "type is ", args.type
        find_alt_r1 = glob.glob(args.pathway_r1 + '/' +  '*qnormed_*' + str(args.threshold_pre) + '*' + '_' + args.chr + '*below0.5.npz')[0]
        find_alt_r2 = glob.glob(args.pathway_r2 + '/' +  '*qnormed_*' + str(args.threshold_pre) + '*' + '_' +args.chr + '*below0.5.npz')[0]
        find_alt_r3 = glob.glob(args.pathway_r3 + '/' +  '*qnormed_*' + str(args.threshold_pre) + '*' + '_' + args.chr + '*below0.5.npz')[0]
        find_alt_r4 = glob.glob(args.pathway_r4 + '/' +  '*qnormed_*' + str(args.threshold_pre) + '*' + '_' + args.chr + '*below0.5.npz')[0]
        find_alt_r5 = glob.glob(args.pathway_r5 + '/' +  '*qnormed_*' + str(args.threshold_pre) + '*' + '_' + args.chr + '*below0.5.npz')[0]
        print "grabbing file: ",find_alt_r1
        print "grabbing file: ",find_alt_r2
        print "grabbing file: ",find_alt_r3
        print "grabbing file: ",find_alt_r4
        print "grabbing file: ",find_alt_r5

        matrix_rep1 = scipy.sparse.load_npz(find_alt_r1)
        matrix_rep1 = matrix_rep1.tocsr()
        matrix_rep2 = scipy.sparse.load_npz(find_alt_r2)
        matrix_rep2 = matrix_rep2.tocsr()
        matrix_rep3 = scipy.sparse.load_npz(find_alt_r3)
        matrix_rep3 = matrix_rep3.tocsr()
        matrix_rep4 = scipy.sparse.load_npz(find_alt_r4)
        matrix_rep4 = matrix_rep4.tocsr()
        matrix_rep5 = scipy.sparse.load_npz(find_alt_r5)
        matrix_rep5 = matrix_rep5.tocsr()
       
    matrix_rep1.data = matrix_rep1.data[np.isnan(matrix_rep1.data) == False]
    matrix_rep1.data = matrix_rep1.data[matrix_rep1.data > args.threshold]
    matrix_rep2.data = matrix_rep2.data[np.isnan(matrix_rep2.data) == False]
    matrix_rep2.data = matrix_rep2.data[matrix_rep2.data > args.threshold]
    matrix_rep3.data = matrix_rep3.data[np.isnan(matrix_rep3.data) == False]
    matrix_rep3.data = matrix_rep3.data[matrix_rep3.data > args.threshold]
    matrix_rep4.data = matrix_rep4.data[np.isnan(matrix_rep4.data) == False]
    matrix_rep4.data = matrix_rep4.data[matrix_rep4.data > args.threshold]
    matrix_rep5.data = matrix_rep5.data[np.isnan(matrix_rep5.data) == False]
    matrix_rep5.data = matrix_rep5.data[matrix_rep5.data > args.threshold]


    matrix_rep1.data = matrix_rep1.data + 1
    matrix_rep2.data = matrix_rep2.data + 1
    matrix_rep3.data = matrix_rep3.data + 1
    matrix_rep4.data = matrix_rep4.data + 1
    matrix_rep5.data = matrix_rep5.data + 1

    final_matrix_rep1 = np.log(matrix_rep1.data)
    final_matrix_rep2 = np.log(matrix_rep2.data)
    final_matrix_rep3 = np.log(matrix_rep3.data)
    final_matrix_rep4 = np.log(matrix_rep4.data)
    final_matrix_rep5 = np.log(matrix_rep5.data)

    #bins = np.linspace(0, final_matrix_rep1.max(), 100)
    bins = np.linspace(0, 4, 100)

    

    pyplot.hist(final_matrix_rep1, bins, normed = True,alpha=0.25, label='merge 0min')
    pyplot.hist(final_matrix_rep2, bins, normed = True,alpha=0.25, label='merge 25min')
    pyplot.hist(final_matrix_rep3, bins, normed = True,alpha=0.25, label='merge 60min')
    pyplot.hist(final_matrix_rep4, bins, normed = True,alpha=0.25, label='merge 120min')
    pyplot.hist(final_matrix_rep5, bins, normed = True,alpha=0.25, label='merge 240min')

    pyplot.legend(loc='upper right')
    #xscale(scale)
   # xtics(scale)

    pyplot.title('Normalized histogram of ' + args.chr + ' after ' + args.type)
    pyplot.xlabel('logged counts')
    if args.type == 'kr' or args.type == 'iced':
        #command_remove = 'rm ' + args.pathway_r1 + '/' + args.chr + '_' + args.type + '_' + args.exp + '_' + str(args.resolution) + '_rep1_rep2_lowest_histogram.png'
        #os.system(command_remove)
        pyplot.savefig(args.pathway_r1 + '/' + args.chr + '_' + str(args.threshold) + '_' + args.type + '_' + str(args.resolution) + '_allmerge_lowest_histogram.png')
    else:
        #command_remove = 'rm ' + args.pathway_r1 + '/' + args.chr + '_' + str(args.threshold) + '_' + args.type + '_' + args.exp + '_' + str(args.resolution) + '_rep1_rep2_lowest_histogram.png'
        #os.system(command_remove)
        pyplot.savefig(args.pathway_r1 + '/' + args.chr + '_' + str(args.threshold) + '_' + args.type + '_' + str(args.resolution) + '_allmerge_lowest_histogram_below0.5.png')

if __name__ == "__main__":
    main()
    
