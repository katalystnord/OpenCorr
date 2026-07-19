/*
 Local build-verification test for NelderMead, issue #17 (Wave 2, stagnant-
 convergence false positive).

 The old "stagnant" guard declared convergence the moment rtol (the relative
 spread of cost across the simplex) happened to land within `tiny` of the
 PREVIOUS iteration's rtol -- a single coincidental near-match, not evidence
 of an actual stall. rtol only measures spread WITHIN the current simplex,
 which can legitimately be small on a very shallow slope that's still far
 (in absolute parameter terms) from the true minimum -- exactly the shape of
 the objective below: a shallow, uniformly-sloped paraboloid whose true
 minimum sits well outside typical starting points and simplex sizes, so
 progress toward it is real but slow, and rtol drifts slowly enough that two
 consecutive samples landing within `tiny` of each other by pure luck is a
 realistic outcome, not a contrived one. The fix requires a run of several
 consecutive non-improving iterations before agreeing the search has
 actually stalled, not just one lucky match.

 Found empirically: a random search over starting points/deltas against a
 standalone reimplementation of both the old and new stagnation logic
 (single-match vs. patience-based) turned up this exact case diverging
 sharply -- old declares convergence after 2 iterations, 141 units from the
 true minimum; new reaches it exactly.
*/

#include <iostream>
#include <cmath>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	//shallow paraboloid, true minimum at (100, 100), f=0
	auto shallow_paraboloid = [](const std::vector<float>& v) -> float
	{
		float x = v[0] - 100.f, y = v[1] - 100.f;
		return 0.0001f * (x * x + y * y);
	};

	cout << "=== Shallow paraboloid: rtol can plateau far from the true minimum ===" << endl;
	cout << "  minimum at (100, 100), f=0" << endl;

	std::vector<float> variables = { 1.59575987f, -4.02070427f }; //found via random search to trigger the false positive
	std::vector<float> deltas = { 0.909981668f, 0.860424221f };
	int iterations_used = 0;
	float final_cost = 0.f;

	NelderMead nm(2000, 1e-8f); //tight tolerance -- the only realistic way to stop early is the stagnation guard, not rtol<tolerance
	bool converged = nm.minimize(variables, deltas, shallow_paraboloid, iterations_used, final_cost);

	cout << "  result: x=" << variables[0] << " y=" << variables[1]
		<< " final_cost=" << final_cost << " iterations=" << iterations_used
		<< " converged=" << (converged ? "true" : "false") << endl;

	float error = std::sqrt((variables[0] - 100.f) * (variables[0] - 100.f) + (variables[1] - 100.f) * (variables[1] - 100.f));
	cout << "  distance from true minimum: " << error << endl;

	//the actual regression: convergence must not be declared while still far from the
	//minimum -- a real (not coincidental) stall. Pre-fix, this case stops after 2
	//iterations at distance ~141; post-fix it reaches the true minimum exactly.
	bool not_premature = error < 1.f;
	cout << "  " << (not_premature ? "PASS" : "FAIL")
		<< ": did not declare convergence prematurely on a coincidental rtol plateau" << endl;
	if (!not_premature) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
