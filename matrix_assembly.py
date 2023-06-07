import argparse
import os
import pandas as pd
import numpy as np
import scipy.sparse as sparse
from iced.normalization import ICE_normalization

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('chrom')
    parser.add_argument('resolution', type=int)
    #parser.add_argument('nz_threshold', type=int)
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    parser.add_argument('-o', '--outdir', default='.')
    parser.add_argument('-c', '--chunksize', default=1000000, type=int)
    parser.add_argument('-n', '--num_chunks', type=int)
    args = parser.parse_args()
    chrom = args.chrom
    resolution = args.resolution
    filenames = args.filenames
    outdir = args.outdir
    chunksize = args.chunksize
    num_chunks = args.num_chunks
    #nz_threshold = args.nz_threshold

    rows = []
    cols = []
    vals = []

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for filename in filenames:
        reader = pd.read_csv(
            filename, sep='\t', chunksize=chunksize, usecols=[1, 2, 4, 5],
            names=['chrom1', 'start1', 'chrom2', 'start2'])

        for idx, chunk in enumerate(reader):
            if num_chunks is not None and idx >= num_chunks:
                break

            print('processing chunk number %s' % idx)

            chunk['start1'] = (chunk['start1'] / resolution).astype(int)
            chunk['start2'] = (chunk['start2'] / resolution).astype(int)
            chunk.drop(chunk[chunk['chrom1'] != chrom].index, inplace=True)
            chunk.drop(chunk[chunk['chrom2'] != chrom].index, inplace=True)
            rows.append(chunk['start1'])
            cols.append(chunk['start2'])
            vals.append(np.ones(len(chunk.index)))


    row = np.concatenate(rows)
    col = np.concatenate(cols)
    val = np.concatenate(vals)

    size = max(np.max(row), np.max(col)) + 1
    matrix = sparse.coo_matrix((val, (row, col)), shape=(size, size))

    #sparse.save_npz('%s/%s_%d_raw.npz' % (outdir,chrom,resolution), matrix)
    sparse.save_npz('%s/%s_raw.npz' % (outdir,chrom), matrix)

    matrix = matrix.tocsr()
    bad_rows = []
    
    ######
    #for i in range(0,len(matrix.indptr)-1): #check number of non-zeros
    #    if matrix.indptr[i+1] - matrix.indptr[i] <= nz_threshold:
    #        bad_rows.append(i)

    #for i in range(0,len(bad_rows)):
    #    matrix.data[matrix.indptr[bad_rows[i]]:matrix.indptr[bad_rows[i]+1]] = 0

    #sparse.save_npz('%s/%s_%d_%d_after_raw.npz' % (outdir,chrom,resolution,nz_threshold), matrix)
    #######

    #matrix = matrix.tocoo()

    #NOTE: turned off ice stage at 60min - not used and eats up CPU time, KR script instead
    #balanced, bias = ICE_normalization(matrix, output_bias=True, verbose=2)
    #sparse.save_npz('%s/%s_%d_iced.npz' % (outdir, chrom,resolution), balanced)
    #np.savetxt('%s/%s_%d_bias.txt' % (outdir, chrom,resolution), bias)


if __name__ == '__main__':

    main()
