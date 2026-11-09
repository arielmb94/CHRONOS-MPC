# CHRONOS documentation development checklist

## Scope

Develop the documentation website around exactly three groups: **Home**,
**User Guide**, and **API Reference**. This roadmap lives at the repository
root, outside the published `docs/` directory.

Existing reference pages should be reviewed and connected to
the new guide, rather than recreated.

## Writing and navigation rules

### Writing

- Assume readers understand MPC and MATLAB but are unfamiliar with CHRONOS
  and LPV-MPC.
- Keep explanations brief and concrete. Prefer short paragraphs, equations,
  tables, and focused code excerpts.
- Explain CHRONOS-specific behavior without repeating general MPC theory.
- Introduce notation and dimensions before use.
- Distinguish exact representations from approximations.
- Follow canonical contracts and source API help; flag unresolved
  discrepancies rather than inventing behavior.

### Independent pages

- Begin each chapter with a brief objective, prerequisites, and expected
  outcome.
- Define the notation needed locally and link directly to any prerequisite
  explanation.
- Avoid requiring readers to complete previous chapters.
- Keep each authoritative explanation in one place; use brief reminders and
  links elsewhere.
- Link to API reference pages instead of duplicating complete function
  documentation.
- Use only one complete worked application initially: a basic two-tank
  controller.

### Navigation and organization

- Every dedicated page must belong to one of the three groups below.
- Use one Markdown file per dedicated chapter. The listed sections become
  headings within that file.
- The sidebar must show **Home**, **User Guide**, and **API Reference**, in
  that order.
- Under User Guide, list the chapters in the order specified below.
- Begin the nonlinear-to-LPV chapter with a short explanation of the target
  CHRONOS formulation before deriving model representations.
- Explain the purpose and timing of coefficient changes in the lifecycle
  chapter; explain array layouts, stage mappings, and update calls in the
  online chapter. Cross-link these explanations.
- Under API Reference, retain **Offline Initialization** and **Online API**,
  with their existing function pages beneath them.
- Give each content page with multiple sections a linked table of contents
  near the beginning.
- Use descriptive titles, logical heading levels, and stable section anchors.
- Provide previous/next navigation and contextual links between related
  pages.
- Keep group index pages brief and ensure every published page is reachable.
- Link from Home directly to the worked tutorial, even though it appears
  later in the User Guide.

## Grouped development checklist

Each row identifies a dedicated page or existing page family and its
sections. Mark its box complete only after the content and navigation have
been reviewed.

### Group 1 — Home

| Dedicated page | Sections |
| --- | --- |
| ☐ Home | Capabilities and limitations; requirements; installation; controller workflow; links to the User Guide, first controller, and API Reference |

### Group 2 — User Guide

| Dedicated page | Sections |
| --- | --- |
| ☐ User Guide index | Ordered chapter links; direct link to the complete first-controller tutorial |
| ☐ From nonlinear problems to LPV-MPC | Brief introduction to the target CHRONOS formulation; signals, notation, and dimensions; nonlinear dynamics and signal equations; scheduling variables; exact representations and local approximations; affine terms and operating-point coordinates; discretization; mapping objectives and constraints to the supported formulation |
| ☐ The CHRONOS controller lifecycle | Define, build, update, solve, apply; offline and online responsibilities; fixed structure versus changing numerical data; why and when coefficients change; changes between samples versus variation across a prediction horizon; connecting scheduling variables to model updates |
| ☐ Defining and building your controller | Dynamics and outputs; known disturbances; tracking, control-effort, and control-increment costs; terminal costs; custom costs and constraints; hard and soft constraints; slack penalties; building the controller |
| ☐ Running and updating your controller | Physical stages and MATLAB indexing; initial and terminal stages; vector columns and matrix pages; constant, partial, and complete horizon data; current state and actually applied previous input; references and known inputs; model, cost, and constraint update functions; update restrictions; solving and applying the first control; retaining solver state; predictions and solver outputs |
| ☐ Putting it all together: your first controller | Basic two-tank problem; model and dimensions; controller definition; building; closed-loop execution; expected results |
| ☐ Controller tuning and troubleshooting | Sample time and horizon; signal scaling; cost weights; slack penalties; solver settings and stopping; tracking problems; constraint violations; ineffective updates; convergence problems |

### Group 3 — API Reference

Existing reference pages will be reviewed and connected to the new guide,
rather than recreated.

| Dedicated page or page family | Sections or review tasks |
| --- | --- |
| ☐ API Reference index | Links to Offline Initialization and Online API |
| ☐ Offline Initialization index and existing function pages | Initialization and building functions; complete function coverage; source-derived help; links to relevant User Guide sections |
| ☐ Online API index and existing function pages | `mpc_solve` first, followed by update functions; complete function coverage; source-derived help; links to relevant User Guide sections |

## Static completion checks

- [ ] Confirm that every checklist entry has an explicit group and navigation
  position.
- [ ] Confirm coverage of nonlinear-to-LPV translation, API mapping, time
  variation, online operation, and tuning.
- [ ] Confirm the chapter order: mathematical representation, controller
  lifecycle, offline definition, online operation, worked controller, tuning.
- [ ] Confirm that the nonlinear-to-LPV chapter introduces the target CHRONOS
  formulation before deriving model representations.
- [ ] Confirm that the lifecycle chapter explains why and when coefficients
  change, while the online chapter explains array layouts, stage mappings,
  and update calls, with cross-links between the explanations.
- [ ] Confirm that the tutorial follows the online-operation chapter and
  precedes tuning.
- [ ] Preserve existing API documentation requirements.

For this roadmap-only change, do not run Sphinx, MATLAB, Octave, examples,
or tests; do not create commits or publish changes.
