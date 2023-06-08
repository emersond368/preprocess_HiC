import pandas as pd
import numpy as np
import scipy.stats as stats
import multiprocessing

def qnorm(data,tie='min'):

    #not nan safe
    data[np.isnan(data)] = 0
    ranks = [stats.rankdata(data[:, i], method=tie)
             for i in range(data.shape[1])]
    workspace = np.array([np.sort(data[:, i])
                          for i in range(data.shape[1])], dtype=float).T
    averages = np.mean(workspace, axis=1)
    for i in range(data.shape[1]):
        workspace[:, i] = averages[ranks[i] - 1]
    return workspace

def qnorm_reducesort(data,tie='min'):

    #not nan safe
    data[np.isnan(data)] = 0
    ranks = [stats.rankdata(data[:, i], method=tie)
             for i in range(data.shape[1])]

    workspace = np.zeros(data.shape)

    #initial fill of sorted values
    for i in range(data.shape[1]):
        workspace[ranks[i]-1,i] = data[:,i]
    
    #fill missing ties left as zero with pandas fillfoward
    df = pd.DataFrame(workspace)
    for i in range(data.shape[1]):
        df[i] = df[i].replace(to_replace=0, method='ffill')
    
    workspace = np.array(df)

    averages = np.mean(workspace, axis=1)
    for i in range(data.shape[1]):
        workspace[:, i] = averages[ranks[i] - 1]

    return workspace

def qnorm_averages(data):
    workspace = np.array([np.sort(data[:, i])
                        for i in range(data.shape[1])], dtype=float).T
    averages = np.mean(workspace, axis=1)

    return averages

def qnorm_averages_parrallel(data):
    "assumes 10 condition data"

    input = [data[:,0:2],data[:,2:4],data[:,4:6],data[:,6:8],data[:,8:10]]  
    pool = multiprocessing.Pool(20)
    results = pool.map(qnorm_averages,input)
  
    overall_average = np.average(results)

    print "overall_average.shape = ",overall_average.shape

    return overall_average

def main():
    assert np.allclose(
        qnorm(np.array([[5, 4, 3],
                        [2, 1, 4],
                        [3, 4, 6],
                        [4, 2, 8]])),
        np.array([[5.66666667, 4.66666667, 2.        ],
                  [2.        , 2.        , 3.        ],
                  [3.        , 4.66666667, 4.66666667],
                  [4.66666667, 3.        , 5.66666667]])
    )
    assert np.allclose(
        qnorm(np.array([[5, np.nan,      3],
                        [2,      1,      4],
                        [np.nan, 4,      6],
                        [4,      2, np.nan]])),
        np.array([[5.        ,     np.nan, 2.        ],
                  [2.        , 2.        , 3.33333333],
                  [    np.nan, 5.        , 5.        ],
                  [3.33333333, 3.33333333,     np.nan]]),
        equal_nan=True
    )


if __name__ == '__main__':
    main()
