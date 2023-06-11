import h5py
import multiprocessing as mp
import cooler
import scipy.sparse
import argparse
import numpy as np

def main():

    parser = argparse.ArgumentParser()    
    parser.add_argument('cooler_input', type=str, help="cooler file")
    parser.add_argument('chr', type=str, help="chromosome")
    parser.add_argument('resolution', type=str, help="resolution")
    args = parser.parse_args()

    #old cooler files (ES and HFF) were powers of 2
    c = cooler.Cooler(args.cooler_input + '::resolutions/' + args.resolution)
    f = h5py.File(args.cooler_input, 'r')

    resolution = c.info['bin-size']
    print("resolution = ",resolution)
    print(args.chr)

    balanced = c.matrix(balance=True).fetch(args.chr)*f['resolutions'][args.resolution]['bins']['weight'].attrs['scale']
    raw = c.matrix(balance=False).fetch(args.chr)
    print("length at balanced = ", balanced.shape[0])
    Sbalanced = scipy.sparse.csr_matrix(balanced)
    print("length at Sbalanced = ", Sbalanced.shape[0])
    Sraw = scipy.sparse.csr_matrix(raw)

    #print out the balanced and raw matrices as sparse .npz matrices
    if "mcool" in args.cooler_input:
        scipy.sparse.save_npz(args.cooler_input[:-6] + '_balanced_' + args.chr + '_'+ args.resolution + '.npz',Sbalanced)
        scipy.sparse.save_npz(args.cooler_input[:-6] + '_raw_' + args.chr + '_'+ args.resolution + '.npz',Sraw)
    else: #assume .cool in name
        scipy.sparse.save_npz(args.cooler_input[:-5] + '_balanced_' + args.chr + '_'+ args.resolution + '.npz',Sbalanced)
        scipy.sparse.save_npz(args.cooler_input[:-5] + '_raw_' + args.chr + '_'+ args.resolution + '.npz',Sraw)

    df = c.bins()[:]


    #print out the bias vector for converting between balanced and raw values
    vector= 1/(df[df["chrom"] == args.chr]["weight"].values*np.sqrt(f['resolutions'][args.resolution]['bins']['weight'].attrs['scale']))

    if "mcool" in args.cooler_input:
        np.savetxt(args.cooler_input[:-6] + '_bias_' + args.chr + '_'+ args.resolution + '.txt',vector)
    else:  #assume .cool in name
        np.savetxt(args.cooler_input[:-5] + '_bias_' + args.chr + '_'+ args.resolution + '.txt',vector)

    #QC

    #checking direction of conversion between raw and balanced (i.e. balanced*bias_1*bias_2 = raw in the cooler files analyzed)
    #if args.chr == 'chr16':
        #balanced_instance = balanced[2000,2100]
        #bias_1 = vector[2000]
        #bias_2 = vector[2100] 
        #raw_instance = raw[2000,2100]

        #print("balanced_instance: ",balanced_instance)
        #print("bias_1: ",bias_1)
        #print("bias_2: ",bias_2)
        #print("raw_instance: ",raw_instance)

        ##assert(raw_instance == np.round(balanced_instance*bias_1*bias_2,6))

if __name__ == '__main__':
    main()
    
    
