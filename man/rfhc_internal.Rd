\name{rfhc_internal}

\alias{rfhc_sample_GMRF}
\alias{rfhc_sample_conditional_GMRF}
\alias{rfhc_blockGMRF}
\alias{rfhc_lMHData}
\alias{rfhc_lPL_WithPL}
\alias{rfhc_interaction_HC}
\alias{rfhc_interaction_HS}
\alias{rfhc_rlambert}
\alias{rfhc_makeQ}
\alias{rfhc_makeQ.full}
\alias{rfhc_makeQtrue}
\alias{rfhc_SOLVE_LAB}
\alias{rfhc_SOLVE_LtAB}
\alias{rfhc_SOLVE_Lvb}
\alias{rfhc_SOLVE_Ltvb}
\alias{rfhc_SOLVE_QAB}
\alias{rfhc_Q2L}
\alias{rfhc_SOLVE_INVQ}
\alias{rfhc_LOGDET}
\alias{rfhc_insideoutside}
\alias{rfhc_node2grid}
\alias{rfhc_grid2node}
\alias{rfhc_xy2grid}
\alias{rfhc_xy2node}
\alias{rfhc.vgram}
\alias{rfhc_modify_pp}
\alias{rfhc_fit_to_im}




\title{Internal Functions}
\description{
This is a list of internally used functions of package \code{\link{rfhc}}.

List provided for transparency.
}
\details{
Simulation:\cr
----------

rfhc_sample_GMRF: used by e.g. simulateGMRF\cr
rfhc_sample_conditional_GMRF: used by e.g. simulate_condGMRF


fit:\cr
----

rfhc_blockGMRF: block update function,for hyper parameters and prediction\cr
rfhc_lMHData\cr
rfhc_lPL_NoPL\cr
rfhc_lPL_WithPL\cr
rfhc_interaction_HC\cr
rfhc_interaction_HS\cr
rfhc_rlambert: Sample Lambert-distribution


matrices:\cr
---------

rfhc_makeQ\cr
rfhc_makeQtrue\cr
rfhc_SOLVE*: solve linear equations using Cholesky factorisation\


grid:\cr
-----

rfhc_insideoutside: make the extended grid and the node representation\cr
rfhc_node2grid: swap between columnwise vector and 2d-grid representation\cr
rfhc_grid2node\cr
rfhc_xy2grid\cr
rfhc_xy2node


other:\cr
-----

rfhc.vgram: estimate variogram using \code{\link{fields}} package.\cr
rfhc_modify_pp: Add z-dimension and rename marks etc, for compatibility with c-functions.\cr
rfhc_fit_to_im

}
\seealso{
\code{\link{rfhc}}.
}
\keyword{rfhc}