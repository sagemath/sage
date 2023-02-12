scip_sdp: Mixed integer semidefinite programming plugin for SCIP
================================================================

Description
-----------

SCIP-SDP allows to solve MISDPs using a nonlinear branch-and-bound
approach or a linear programming cutting-plane approach.

- In the first case (the default), the semidefinite programming (SDP)
  relaxations are solve using interior-point SDP-solvers.

- In the second case, cutting planes based on eigenvector are
  generated.

SCIP-SDP is based on the branch-and-cut framework SCIP. In addition to
providing a constraint handler for SDP-constraints and a relaxator to
solve continuous SDP-relaxations using interior-point solvers,
SCIP-SDP adds several heuristics and propagators to SCIP.

License
-------

Apache 2.0


Upstream Contact
----------------

http://www.opt.tu-darmstadt.de/scipsdp/

https://github.com/scipopt/SCIP-SDP
