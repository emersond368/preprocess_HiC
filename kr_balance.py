#import sys
#sys.path.insert(0, '/project/jcreminslab/tgg_projects_2/bonev')

import argparse
import subprocess
import os

import pandas as pd
import scipy.sparse as sparse

from config import KR_JAR


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('chrom')
    parser.add_argument('nz_threshold',type = int)
    parser.add_argument('-r', '--resolution', type=int, default=10000)
    args = parser.parse_args()

    directory = args.directory
    chrom = args.chrom
    resolution = args.resolution

    if args.resolution == 10000:
        raw_npz_file = '%s/%s_raw.npz' % (directory, chrom)
        raw_txt_file = '%s/%s_raw.txt' % (directory, chrom)
        raw_txt_temp_file = '%s/%s_raw_temp.txt' % (directory, chrom)
        hic_file = '%s/%s.hic' % (directory, chrom)
        bias_file = '%s/%s_kr.bias' % (directory, chrom)
        kr_txt_file = '%s/%s_kr.txt' % (directory, chrom)
        kr_npz_file = '%s/%s_kr.npz' % (directory, chrom)
        old_hic_file = '%s/%s_%d_kr.hic' % (directory, chrom, resolution)
        old_bias_file = '%s/%s_%d_kr.bias' % (directory, chrom, resolution)
        old_kr_txt_file = '%s/%s_%d_kr.txt' % (directory, chrom, resolution)
        old_kr_npz_file = '%s/%s_%d_kr.npz' % (directory, chrom, resolution)
        old_raw_txt_file = '%s/%s_%d_raw.txt' % (directory, chrom, resolution)
        
        if os.path.isfile(raw_txt_file):
            os.remove(raw_txt_file)
        if os.path.isfile(kr_npz_file):
            os.remove(kr_npz_file)
        if os.path.isfile(hic_file):
            os.remove(hic_file)
        if os.path.isfile(old_hic_file):
            os.remove(old_hic_file)
        if os.path.isfile(old_bias_file):
            os.remove(old_bias_file)
        if os.path.isfile(old_kr_txt_file):
            os.remove(old_kr_txt_file)
        if os.path.isfile(old_kr_npz_file):
            os.remove(old_kr_npz_file)
        if os.path.isfile(old_kr_txt_file):
            os.remove(old_kr_txt_file)

    else:
        raw_npz_file = '%s/%s_%d_raw.npz' % (directory, chrom,resolution)
        raw_txt_file = '%s/%s_%d_raw.txt' % (directory, chrom,resolution)
        raw_txt_temp_file = '%s/%s_%d_raw_temp.txt' % (directory, chrom,resolution)
        hic_file = '%s/%s_%d.hic' % (directory, chrom,resolution)
        bias_file = '%s/%s_%d_kr.bias' % (directory, chrom,resolution)
        kr_txt_file = '%s/%s_%d_kr.txt' % (directory, chrom,resolution)
        kr_npz_file = '%s/%s_%d_kr.npz' % (directory, chrom,resolution)
        old_hic_file = '%s/%s.hic' % (directory, chrom)
        old_bias_file = '%s/%s_kr.bias' % (directory, chrom)
        old_kr_txt_file = '%s/%s_kr.txt' % (directory, chrom)
        old_kr_npz_file = '%s/%s_kr.npz' % (directory, chrom)
        old_raw_txt_file = '%s/%s_raw.txt' % (directory, chrom)

        if os.path.isfile(old_hic_file):
            os.remove(old_hic_file)
        if os.path.isfile(old_bias_file):
            os.remove(old_bias_file)
        if os.path.isfile(old_kr_txt_file):
            os.remove(old_kr_txt_file)
        if os.path.isfile(old_kr_npz_file):
            os.remove(old_kr_npz_file)
        if os.path.isfile(raw_txt_file):
            os.remove(raw_txt_file)
        if os.path.isfile(hic_file):
            os.remove(hic_file)
        if os.path.isfile(kr_npz_file):
            os.remove(kr_npz_file)
        if os.path.isfile(bias_file):
            os.remove(bias_file)
        if os.path.isfile(kr_txt_file):
            os.remove(kr_txt_file)
        if os.path.isfile(old_raw_txt_file):
            os.remove(old_raw_txt_file)

    matrix = sparse.load_npz(raw_npz_file).tocsr()
    size = matrix.shape[0]
    indices = matrix.indices
    indptr = matrix.indptr
    data = matrix.data

    ##turn on row removal
    bad_rows = []
    
    for i in range(0,len(matrix.indptr)-1): #check number of non-zeros
        if matrix.indptr[i+1] - matrix.indptr[i] <= args.nz_threshold:
            bad_rows.append(i)

    print "number of bad rows < " +  str(args.nz_threshold) + ": ", len(bad_rows)
    for i in range(0,len(bad_rows)):
        matrix.data[matrix.indptr[bad_rows[i]]:matrix.indptr[bad_rows[i]+1]] = 0 

    
    for i in range(size):
        if i % 100 == 0:
            print('writing row %i/%i' % (i, size))
        df = pd.DataFrame({'i': i, 'j': indices[indptr[i]:indptr[i+1]],
                           'v': data[indptr[i]:indptr[i+1]]})
        df['str1'] = 0
        df['chr1'] = chrom
        df['pos1'] = df['i'] * resolution
        df['frag1'] = 0
        df['str2'] = 1
        df['chr2'] = chrom
        df['pos2'] = df['j'] * resolution
        df['frag2'] = 1
        df['score'] = df['v']
        df.drop(['i', 'j', 'v'], axis=1, inplace=True)
        df.to_csv(raw_txt_file, sep='\t', header=False, index=False, mode='a')
    
    cmd = 'sort -u %s > %s'%(raw_txt_file,raw_txt_temp_file)
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'mv %s %s'%(raw_txt_temp_file,raw_txt_file)
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = '%s pre -d -c %s -r %i %s %s mm9' % \
        (KR_JAR, chrom, resolution, raw_txt_file, hic_file)
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = '%s dump observed KR %s %s %s BP %i %s' % \
        (KR_JAR, hic_file, chrom, chrom, resolution, kr_txt_file)
    print(cmd)
    subprocess.call(cmd, shell=True)

    print('saving %s' % kr_npz_file)
    df = pd.read_csv(kr_txt_file, sep='\t', names=['i', 'j', 'v'])
    df['i'] = (df['i'] / resolution).astype(int)
    df['j'] = (df['j'] / resolution).astype(int)
    matrix_kr = sparse.coo_matrix((df['v'], (df['i'], df['j'])),
                                  shape=(size, size))
    sparse.save_npz(kr_npz_file, matrix_kr)

    cmd = '%s dump norm KR %s %s BP %i %s' % \
        (KR_JAR, hic_file, chrom, resolution, bias_file)
    print(cmd)
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
