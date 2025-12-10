# 3-AFC Metacognitive Sensitivity Task  
*A Psychtoolbox-powered random-dot kinematogram for measuring how well confidence tracks visual-decision accuracy.*

---

## Overview  
This repository delivers a complete MATLAB experiment that replicates and extends Zehetleitner & Rausch (2013).  
Participants view a circular aperture of moving dots, decide whether the dominant motion is left, right, or absent, and then rate their confidence.  
Pupil-size changes can be recorded in parallel via any eye-tracker that streams to a serial port (optional).  
The code exports trial-wise data ready for Type-2 signal-detection and meta-d′ analysis.

---

## What you get  
- **9 coherence levels** (0 % – 100 %) fully randomised across 441 trials  
- **3-AFC response**: left-arrow / right-arrow / “N” for pure noise  
- **6-point Likert confidence scale** (mouse-driven, 1 = very low → 6 = very high)  
- **Adaptive Gabor mask** (45°, 1-s fade) hides stimulus onset  
- **Automatic breaks** every block; run-time ≈ 55 min  
- **Instant CSV output** (`PXXX_Data.csv`) with fields:  
  `participantNumber, trial, block, coherence, direction, response, correct, confidence, RT, age, sex`

---

## Quick start  
1. Clone or download the repo  
2. Open MATLAB, `cd` into the folder  
3. Run `Metacognition_Test.m`  
4. Enter participant ID, age, sex when prompted  
5. Follow on-screen instructions; data file appears in the same folder once the task ends

---

## Requirements  
- MATLAB R2020b or later  
- Psychtoolbox-3 (free)  
- Windows / macOS / Linux with ≥ 120 Hz display recommended  
- (Optional) eye-tracker streaming pupil diameter to a COM port

---

## File map  
| File | Purpose |
|------|---------|
| `Metacognition_Test.m` | Main experiment driver |
| `movingDots.m` | Stimulus routine (coherent + random dots, Gabor overlay) |
| `getConfidenceRating.m` | Interactive Likert scale (mouse hover + click) |
| `generateCircularGabor.m` | Creates oriented Gabor patch texture |
| `constructContingencyTables.m` | Builds Type-2 SDT table and writes CSV |
| `angle2pix.m` | Converts visual angle → screen pixels |
| `displayInstructions.m` | Full-screen instruction slides |

---

## Analysis tips  
- **Type-2 hits / misses / false-alarms / correct-rejections** are already coded by `constructContingencyTables`  
- Feed the CSV into [meta-d′ MATLAB toolbox](https://github.com/metacog/meta-d) or compute classic `Aroc` / `Broc` directly  
- Coherence × confidence interaction? A simple mixed ANOVA (confidence as DV, coherence as within factor) usually reveals the expected monotonic relationship when the task works

---

## Citation  
If you use this code, please cite:  
&gt; Choudra, A. (2024). *Exploring Metacognitive Sensitivity through a 3-AFC Visual Motion Discrimination Task*. MSc Dissertation, University of Reading.

and the original paper:  
&gt; Zehetleitner, M., & Rausch, M. (2013). Being confident without seeing: What subjective measures of visual consciousness are about. *Attention, Perception, & Psychophysics*, 75(7), 1406–1426. https://doi.org/10.3758/s13414-013-0505-2
