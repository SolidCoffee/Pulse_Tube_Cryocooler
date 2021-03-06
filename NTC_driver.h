#ifndef NTC_driver_h
#define NTC_driver_h

#include "Arduino.h"

int R_RTD[]={329,
330,
332,
333,
335,
336,
337,
339,
340,
341,
343,
344,
346,
347,
348,
350,
351,
352,
353,
355,
356,
357,
359,
360,
361,
363,
364,
365,
366,
368,
369,
370,
371,
373,
374,
375,
376,
378,
379,
380,
381,
382,
384,
385,
386,
387,
388,
390,
391,
392,
393,
394,
395,
397,
398,
399,
400,
401,
402,
403,
404,
406,
407,
408,
409,
410,
411,
412,
413,
414,
416,
417,
418,
419,
420,
421,
422,
423,
424,
425,
426,
427,
428,
429,
430,
431,
432,
433,
434,
435,
436,
437,
438,
439,
440,
441,
442,
443,
444,
445,
446,
447,
448,
449,
450,
451,
452,
453,
454,
455,
456,
457,
458,
459,
460,
461,
462,
462,
463,
464,
465,
466,
467,
468,
469,
470,
471,
472,
472,
473,
474,
475,
476,
477,
478,
479,
479,
480,
481,
482,
483,
484,
485,
485,
486,
487,
488,
489,
490,
490,
491,
492,
493,
494,
495,
495,
496,
497,
498,
499,
499,
500,
501,
502,
503,
503,
504,
505,
506,
507,
507,
508,
509,
510,
510,
511,
512,
513,
514,
514,
515,
516,
517,
517,
518,
519,
520,
520,
521,
522,
523,
523,
524,
525,
525,
526,
527,
528,
528,
529,
530,
530,
531,
532,
533,
533,
534,
535,
535,
536,
537,
537,
538,
539,
540,
540,
541,
542,
542,
543,
544,
544,
545,
546,
546,
547,
548,
548,
549,
550,
550,
551,
552,
552,
553,
554,
554,
555,
555,
556,
557,
557,
558,
559,
559,
560,
561,
561,
562,
562,
563,
564,
564,
565,
565,
566,
567,
567,
568,
569,
569,
570,
570,
571,
572,
572,
573,
573,
574,
575,
575,
576,
576,
577,
577,
578,
579,
579,
580,
580,
581,
581,
582,
583,
583,
584,
584,
585,
585,
586,
587,
587,
588,
588,
589,
589,
590,
590,
591,
592,
592,
593,
593,
594,
594,
595,
595,
596,
596,
597,
597,
598,
599,
599,
600,
600,
601,
601,
602,
602,
603,
603,
604,
604,
605,
605,
606,
606,
607,
607,
608,
608,
609,
609,
610,
610,
611,
611,
612,
612,
613,
613,
614,
614,
615,
615,
616,
616,
617,
617,
618,
618,
619,
619,
620,
620,
621,
621,
622,
622,
622,
623,
623,
624,
624,
625,
625,
626,
626,
627,
627,
628,
628,
628,
629,
629,
630,
630,
631,
631,
632,
632,
633,
633,
633,
634,
634,
635,
635,
636,
636,
637,
637,
637,
638,
638,
639,
639,
640,
640,
640,
641,
641,
642
};

int T_RTD[]={-152,
-151,
-150,
-149,
-148,
-147,
-147,
-146,
-145,
-144,
-143,
-142,
-142,
-141,
-140,
-139,
-138,
-137,
-136,
-136,
-135,
-134,
-133,
-132,
-131,
-131,
-130,
-129,
-128,
-127,
-126,
-125,
-125,
-124,
-123,
-122,
-121,
-120,
-120,
-119,
-118,
-117,
-116,
-115,
-114,
-114,
-113,
-112,
-111,
-110,
-109,
-109,
-108,
-107,
-106,
-105,
-104,
-103,
-103,
-102,
-101,
-100,
-99,
-98,
-98,
-97,
-96,
-95,
-94,
-93,
-92,
-92,
-91,
-90,
-89,
-88,
-87,
-87,
-86,
-85,
-84,
-83,
-82,
-81,
-81,
-80,
-79,
-78,
-77,
-76,
-75,
-75,
-74,
-73,
-72,
-71,
-70,
-70,
-69,
-68,
-67,
-66,
-65,
-64,
-64,
-63,
-62,
-61,
-60,
-59,
-59,
-58,
-57,
-56,
-55,
-54,
-53,
-53,
-52,
-51,
-50,
-49,
-48,
-48,
-47,
-46,
-45,
-44,
-43,
-42,
-42,
-41,
-40,
-39,
-38,
-37,
-37,
-36,
-35,
-34,
-33,
-32,
-31,
-31,
-30,
-29,
-28,
-27,
-26,
-26,
-25,
-24,
-23,
-22,
-21,
-20,
-20,
-19,
-18,
-17,
-16,
-15,
-15,
-14,
-13,
-12,
-11,
-10,
-9,
-9,
-8,
-7,
-6,
-5,
-4,
-4,
-3,
-2,
-1,
0,
1,
2,
2,
3,
4,
5,
6,
7,
7,
8,
9,
10,
11,
12,
13,
13,
14,
15,
16,
17,
18,
19,
19,
20,
21,
22,
23,
24,
24,
25,
26,
27,
28,
29,
30,
30,
31,
32,
33,
34,
35,
35,
36,
37,
38,
39,
40,
41,
41,
42,
43,
44,
45,
46,
46,
47,
48,
49,
50,
51,
52,
52,
53,
54,
55,
56,
57,
57,
58,
59,
60,
61,
62,
63,
63,
64,
65,
66,
67,
68,
68,
69,
70,
71,
72,
73,
74,
74,
75,
76,
77,
78,
79,
79,
80,
81,
82,
83,
84,
85,
85,
86,
87,
88,
89,
90,
90,
91,
92,
93,
94,
95,
96,
96,
97,
98,
99,
100,
101,
102,
102,
103,
104,
105,
106,
107,
107,
108,
109,
110,
111,
112,
113,
113,
114,
115,
116,
117,
118,
118,
119,
120,
121,
122,
123,
124,
124,
125,
126,
127,
128,
129,
129,
130,
131,
132,
133,
134,
135,
135,
136,
137,
138,
139,
140,
140,
141,
142,
143,
144,
145,
146,
146,
147,
148,
149,
150,
151,
151,
152,
153,
154,
155,
156,
157,
157,
158,
159,
160,
161,
162,
162,
163,
164,
165,
166,
167,
168,
168,
169,
170,
171,
172,
173,
173,
174,
175,
176,
177,
178,
179,
179,
180,
181,
182,
183,
184,
184,
185,
186,
187,
188,
189,
190
};

int i=0;
int Temp_C;
int Temp_F;
int Temp_K;
int RTCval_H;
int RTCval;

#endif
