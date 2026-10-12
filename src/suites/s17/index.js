// Suite 17 — Engineering Economics and Techno-Economic Analysis.
import investmentAppraisal from './investment-appraisal.js';
import capexEstimate from './capex-estimate.js';
import parametricCost from './parametric-cost.js';
import costRisk from './cost-risk.js';
import scheduleRisk from './schedule-risk.js';
import opexLifecycle from './opex-lifecycle.js';
import replacement from './replacement-service-life.js';
import decommissioning from './decommissioning.js';
import productionRevenue from './production-revenue.js';
import renewableLcoe from './renewable-lcoe.js';
import projectFinance from './project-finance.js';
import uncertaintyRisk from './uncertainty-risk.js';
import reliabilityMaintenance from './reliability-maintenance.js';
import optimisationDecision from './optimisation-decision.js';
import realOptions from './real-options.js';
import economicForecasting from './economic-forecasting.js';
import technoEconomic from './techno-economic.js';

export default {
  n: 17,
  modules: [investmentAppraisal, capexEstimate, parametricCost, costRisk, scheduleRisk, opexLifecycle, replacement, decommissioning,
    productionRevenue, renewableLcoe, projectFinance, uncertaintyRisk, reliabilityMaintenance, optimisationDecision, realOptions, economicForecasting, technoEconomic],
};
