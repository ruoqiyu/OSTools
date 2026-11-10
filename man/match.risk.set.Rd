\name{match.risk.set}
\alias{match.risk.set}
\alias{print.os_risk_set_match}
\title{Chronological optimal nonoverlapping matching at treatment decision times}
\usage{
match.risk.set(dat, id, time, z, covariates = character(),
  exact = character(), ncontrol = 1, caliper = Inf, eligible = NULL,
  type = c("initiation", "event"))
\method{print}{os_risk_set_match}(x, ...)
}
\arguments{
\item{x}{An object returned by the corresponding function.}
\item{...}{Additional print arguments (currently unused).}
\item{dat}{Long data, with at most one row per person and decision time.}
\item{id,time,z}{Column names for person ID, ordered numeric or Date decision time, and 0-1 exposure event. In initiation mode z marks initiation, not ongoing treatment status.}
\item{covariates}{Names of pretreatment covariates for rank-Mahalanobis distance. Numeric, logical, factor and character columns are supported. Empty means all permitted distances are zero.}
\item{exact}{Column names requiring exact agreement at the index time.}
\item{ncontrol}{Positive integer number of controls per treated individual.}
\item{caliper}{Maximum square-root rank-Mahalanobis distance. Inf imposes none.}
\item{eligible}{Optional name of a nonmissing logical or 0-1 eligibility column. Eligibility must be defined without consulting later outcomes or later treatment.}
\item{type}{Initiation excludes people already treated before the index time; event compares current event types at each decision, allowing repeated events for unmatched people.}
}
\description{Constructs complete one-treated-to-ncontrol sets in chronological order, using only the current eligible rows to compute distances. Each person appears in at most one matched set.}
\details{
The caller must provide decision-time records and covariates measured
before each decision. The function does not reconstruct histories from
outcomes, interpolate records from other times, or prove that a selected
column is pretreatment. All matching and exact variables must be complete;
outcomes and other unused columns can be missing. In initiation mode the
records must include all prior initiation events relevant to eligibility.
Controls are not selected according to whether they will ever be treated.
A person selected as a control is not reused when later treated, and that
later treated row is reported as already matched.

At each time, rank-Mahalanobis distances are computed from the currently
eligible pool using untied-rank variance adjustment and a generalized
inverse for collinearity. Exact restrictions and the distance caliper are
then applied. Within each time and exact stratum, binary integer assignment
first maximizes the number of complete one-treated-to-ncontrol sets and then
minimizes the total treated-to-control distance at that cardinality. Thus a
treated individual with too few permitted controls is omitted, and competing
calipers or shortages may cause other treated individuals to be omitted.
No partial sets are made. The two optimization stages use Rglpk with sparse
constraints from slam; a solver failure is an error, never a greedy fallback.
Optimality is subject to the solver's numerical tolerances. Input rows are
ordered by ID before assignment, so row permutations do not change the match
in the same solver environment. Equally optimal assignments may differ across
solver versions or platforms.

The optimization is conditional on earlier matches. It does not maximize the
number of sets or minimize distance jointly across all times, and it does not
use future covariates or events to improve a current assignment. Inspect
unmatched.treated and diagnostics and check balance before outcome analysis.

The time column identifies a common risk-set opportunity, such as a
calendar date for initiation or pregnancy number for a birth-type
comparison. When using pregnancy number, include maternal age and relevant
prior history among the covariates/exact variables. This implementation
supports nonoverlapping matched sets; it does not support control reuse or
risk-set sampling weights. Event mode does not itself establish haphazard
assignment of the event type.

Output column names mset, .risk.time and .risk.row are reserved. Original
rows and columns are preserved in data, with those three columns appended.
}
\value{An os_risk_set_match list with data, z (matched exposure indicators),
selected (logical input-row mask), unmatched.treated (row, ID, time, reason),
diagnostics by time, ncontrol, type, algorithm and call. Empty matches
return an empty data frame and diagnostics rather than malformed sets.}
\examples{
d <- expand.grid(id=1:8,time=1:2)
d$z <- as.integer(d$id==2*d$time)
d$history <- d$id %% 2
fit <- match.risk.set(d,"id","time","z",covariates="history")
fit
fit$data
}
\references{
Li, Y. P., Propert, K. J. and Rosenbaum, P. R. (2001).
Balanced risk set matching. Journal of the American Statistical Association
96, 870--882.

Zubizarreta, J. R., Small, D. S. and Rosenbaum, P. R. (2014).
Isolation in the construction of natural experiments.
The Annals of Applied Statistics 8, 2096--2121.
}
