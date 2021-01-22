<a href="http://mit.edu/sparklab/">
  <img align="left" src="http://web.mit.edu/sparklab/assets/images/spark-shard.png" height="80px" alt="sparklab">
</a>

<h1 style="padding-left: 1em; height: 70px;"> Outlier robust estimation </h1>


This repository contains the MATLAB implementation of **GNC (Graduated Non-Convexity)** and **ADAPT (Adaptive Trimming)** described in the following papers:

- Antonante, P., Tzoumas, V., Yang, H., & Carlone, L. (2020). "Outlier-Robust Estimation: Hardness, General-Purpose Algorithms, Experiments, and Guarantees."

```bibtex
@techreport{Antonante20tr-outlierRobustEstimation,
  title = {Outlier-Robust Estimation:
  	Hardness, General-Purpose Algorithms, Experiments, and Guarantees},
  author = {Antonante, P. and Tzoumas, V. and Yang, H. and Carlone, L.},
  hidenote = {in preparation},
  year = {2020}
}
```

- Yang, H., Antonante, P., Tzoumas, V., & Carlone, L. (2020). "Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection". IEEE Robotics and Automation Letters (RA-L), 5(2), 1127–1134.

```bibtex
@article{Yang20ral-GNC,
  author = {Yang, H. and Antonante, P. and Tzoumas, V. and Carlone, L.},
  fullauthor = {Heng Yang, Pasquale Antonante, Vasileios Tzoumas, Luca Carlone},
  title = {Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection},
  volume = {5},
  number = {2},
  pages = {1127--1134},
  pdf = {https://arxiv.org/pdf/1909.08605.pdf},
  journal = {{IEEE} Robotics and Automation Letters ({RA-L})},
  year = {2020}
}
```

- Tzoumas, V., Antonante, P., & Carlone, L. (2019). Outlier-Robust Spatial Perception: Hardness, General-Purpose Algorithms, and Guarantees. IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS).

```bibtex
@inproceedings{Tzoumas19iros-outliers,
  author = {Tzoumas, V. and Antonante, P. and Carlone, L.},
  title = {Outlier-Robust Spatial Perception: Hardness, General-Purpose Algorithms, and Guarantees},
  hidenote = {Extended arxiv version: 1903.11683, \linkToPdf{https://arxiv.org/pdf/1903.11683.pdf}},
  pdf = {https://arxiv.org/pdf/1903.11683.pdf},
  booktitle = {IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS)},
  year = {2019}
}
```

## Quick-start

Open matlab and run

```matlab
setup
```

This will add all the folders to the path.
I suggest to run the script every time you open MATLAB (to prevent path pollution) but you can run `savepath` to save the changes.

Now explore and run the `example.m`.

## GNC (Graduated Non-Convexity) Example

All algorithms provide similar interfaces, let's use GNC as example to solve a **linear regression problem**.
First let's generate a random problem with 100 measurements, 80% of them being outliers:

```matlab
problem = linearRegressionProblem(100, 0.8);
```

Suppose suppose our inlier threshold (epsilon) is set using the chi2 distribution

```matlab
epsilon = chi2inv(0.99, problem.dof)*problem.MeasurementNoiseStd^2
```

We can run GNC simply running

```matlab
[inliers, info] = gnc(problem, @leastSquareNorm2, 'NoiseBound', epsilon);
```

where `@leastSquareNorm2` is the function handle of the non-minimal solver.
The GNC function will return the estimated set of inliers together with other diagnostic information.

## Acknowledgments

This work was partially funded by:

- ARL DCIST CRA W911NF-17-2-0181
- ONR RAIDER N00014-18-1-2828
- Lincoln Laboratory’s “Resilient Perception in Degraded Environments” program.
- Mathworks

## License

[BSD License](LICENSE.BSD)

<div align="center">
  <a href="https://mit.edu">
    <img src="http://mit.edu/sparklab/assets/images/mit.png" width="100" alt="mit">
  </a>
</div>
