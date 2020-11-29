Main program:
huber_boot.m: a Huber regression based inference procedure (Algorithm 1)
huber_panel_boot.m: Huber Robust Multiple Testing procedure (Algorithm 2)

Test program:
test_huber_infer.m
test_panel_huber.m

Misc:
We use minFunc (M. Schmidt. minFunc: unconstrained differentiable multivariate optimization in Matlab. http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.) when inplementing our algorithms


Please cite the the following paper if the package was used

Xi Chen and Wen-Xin Zhou. Robust Inference via Multiplier Bootstrap. Annals of Statistics, 48(3): 1665–1691 2020.

@article{chen2020robust,
  title={Robust Inference via Multiplier Bootstrap},
  author={Xi Chen and Wen-Xin Zhou},
  journal={Annals of Statistics},
  volume={48},
  number={3},	
  year={2020}
}