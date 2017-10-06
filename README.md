# Sliding-window ribbon carving for video condensation

This project is to shorten the length of surveillance video while preserving motion activities in the video. Video ribbon carving is an algorithm based on image seam carving which removes unimportant parts. We adopt a sliding-window approach to process very long video sequences, in principle infinite-length ones.

## Prerequisites

1. MATLAB environment

2. GIT LFS is recommended to download test videos

3. Enough storage to store output videos (Estimation: if storage for current input videos is `X` and max flex parameter is `flexmax`, then `(flexmax + 1) * X` extra space should be enough to store output videos.)


## Running the tests

1. Compile ribboncarvemainC.c in MATLAB:

```
mex ribboncarvemainC.c
```

If it is successfully compiled, it should generate a MATLAB MEX file with same filename.

If it cannot be compiled successfully, you may try rolling back to an earlier commit. The first commit should work fine if it compiles. The original MATLAB code is for reference and could be used for testing as well.

2. Have 'whole_video.avi' and 'whole_cost.avi' files in the same directory. ’whole_video.avi' is the original video and 'whole_cost.avi' is its background subtraction result (cost video). The avi files may require certain format, depending on if aviread() function can decode them properly in MATLAB. You may use sample videos in the repo by downloading via Git LFS.

3. Run scriptC.m in MATLAB.

```
scriptC
```

If it runs successfully, the program will generate output videos (given flexmax = 3): 'flex0_video.avi','flex0_cost.avi',
'flex1_video.avi','flex1_cost.avi',
'flex2_video.avi','flex2_cost.avi',
'flex3_video.avi','flex3_cost.avi'

where 'flex3_video.avi' is our final condensation result.

## Future Improvements

* Use memmove() in carveRibbon() instead of moving pixel values one by one. For doing this we may need to know the offset index in the videos so that we know where we can move a chunk of memory.

* We can skip computing minimum cost of one direction if the other one can proceed. This saves computation time.

* Downsize for speed: we probably don’t have to really downsize the video, but when calculating minimum costs and paths, we can do it every other pixel or so. The downside is that the result may not be the same.

* Carve multiple ribbons at once. For this we may need to redesign path finding algorithm as paths may overlap.

* Parallel processing. For example, calculating vertical and horizontal ribbons in different threads.

* Reduce space footage in dynamic programming.

* Port more code to C when necessary.


## Reference

For more details, please refer to my technical report:
http://iss.bu.edu/data/jkonrad/reports/Wu09-03buece.pdf

## Author

* **Huan-Yu Wu**
