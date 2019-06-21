function dlist=get_dlist()

%function dlist=get_dlist()
%   get dlist: datafile number
%
% see also display_datalist
%
% Corentin Massot
% Cognition and Sensorimotor Integration Lab, Neeraj J. Gandhi
% University of Pittsburgh
% created 11/03/2015 last modified 01/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uday's data VG
%dlist=[1:9]

%'bl_rSCTrack_031115_NoPuff_Overlap.mat'
%20deg

%'bb_rSCTrack_031315_NoPuff_Overlap.mat'
%>15deg

%'bb_lSCTrack_070915_NoPuff_Overlap.mat'
%12deg

%'bb_lSCTrack_071215_NoPuff_Overlap.mat'
%8deg

%'bl_lSCTrack_071415_NoPuff_Overlap.mat'
%11deg

%'bl_lSCTrack_072315_1_NoPuff_Overlap.mat'
%11deg

%'bl_lSCTrack_072315_2_NoPuff_Overlap.mat'
%11deg %issue with VMI

%'bb_lSCTrack_080415_NoPuff_Overlap.mat'
%8deg

%'bb_mSCtrack_121415_NoPuff_Overlap.mat'
%15deg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%022616 2deg ROSTRAL ok sacc tracking
%issue 3 channel 11 lfp noisy
%dlist=9+[3 5]%[1 2 3 4 5 6]

%032116 10deg ok interesting patterns
%dlist=9+[8 9 10 11]%[7 8 9 10 11 12]

%032916 5deg ROSTRAL open RF no sacc tracking
%dlist=9+[14 15 16 17]%[13 14 15 16 17]

%033016 5deg ROSTRAL ok  open RF sacc tracking
%dlist=9+[19 20 21 22]% 21 22];%[18 19 20 21 22]
%dlist=[28]

%%%033116 NOT responding
%%dlist=9+[24 25 26 27]%[23 24 25 26 27]

%062416 issue 28(nch) 29 31 32 because only 1 trial in opposite side!
%%dlist=9+[29]%[28 29 30 31 32]


%120816 10deg 
%dlist=9+34%[33 34]

%120916 >20deg CAUDAL  38 
%dlist=9+[35 36 40];%[35 36 37 38 39 40]

%121216 5deg ROSTRAL ok sacc tracking
%dlist=9+[41 42 44 46];%9+[41 42 43 44 45 46]

%121316 15deg ok
%dlist=9+[47 48 52 54]%[47 48 49 50 51 52 53 54]

%121616 3deg ROSTRAL ok targ latency gradiation
%dlist=[64    65    69    71]%9+[55 56 57 58 59 60 61 62]

%121816 4deg ROSTRAL ok
%dlist=9+[63 64 65 66 67 68]


%020217 10deg / low response / ISSUE Alignement / DISCARDED
%dlist=[78    79    80    81    82    83]

%020817 8deg open RF deep layers maybe tracking only build-up cells
%dlist=[84    85    86    87    88    89    90]
%dlist=88

%031817 7deg / lower responses / LFP look different reversal not at the same place
%dlist=[91 92 93 94 95 96 97]

%032117 7deg / openRf  / all layers very good spreading / small visual burst
%dlist=[98 99 100 102 103 104]%[98 99 100 101 102 103 104]

%032217 7deg /only 3 or 4 channels responding / bad LFP for gap_gapstim and step_basic
%dlist=[105 106 107 109 110 111] %[105 106 107 108 109 110 111]

%032917 7deg / good / 2 depths / ~no visual layers?
%dlist=[112 113 114 115 118 119 120 121 122 123 124]
%dlist=[112 113 114 115 116 117 118 119 120 121 122 123 124]
%dlist=[114]

%033017 7deg / good but some channels were less responsive during expe /
%bad LFP: maybe because of bad cable
%dlist=[125 126 127 128 129 130 131]


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All data / removed discarded data
%dlist=[1:46 48:77 84:131];

% dlist=[ 1     2     3     4     5     6     7     8     9   12      14      20 ...
%         44    50    56    64    72    84    91    98   105   112   115 ...
%     17    18    45    51    57    65    73      86    93   100   107   113   114 ...
%     23:26 28:31]

%dlist=[44    50    56    64    72    84    91    98   105   112   115 ...
%        1     2     3     4     5     6     7     8     9    12    14    20 ...
%       17    18    45    51    57    65    73    86    93   100   107   113   114]

%-56 -105 -9 -57 -107
%dlist=[73 44    50      64    72    84    91    98     112   115 ...
%    1     2     3     4     5     6     7     8       12    14    20 ...
%    17    18    45    51        65      86    93   100      113   114]
      
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all data 
%VMI analysis vg mg step gap

%%%%%%%%%%
%all VG %NOTE add 03/29 03/30?  9 should be removed
%dlist=[  1   2   3     4     5     6     7     8     9   12      14      20 ...
%    44    50    56    64    72    84    91    98   105   112   115 ]


%all VG +28
%dlist=[28  1   2   3     4     5     6     7     8     9   12      14      20 ...
%44    50    56    64    72    84    91    98   105   112   115 ]

%VG no nopuff
%dlist=[  12      14      20 ...
%    44    50    56    64    72    84    91    98   105   112   115 1   2   3     4     5     6     7     8     9]


%dlist=[28 64 86];
%dlist=[64];
%dlist=[86];
%dlist=[64 86];
%dlist=[98 84]
%dlist=84;
%dlist=28;

%RESULTS
%DATA VG -105 -56 -9
%dlist=[ 12      14      20 44    50       64    72    84    91    98     112   115 ...
% 1   2   3     4     5     6     7     8    ]
dlist=[64    72    84    91    98     112   115  12      14      20 44    50        ...
 1   2   3     4     5     6     7     8    ]



% %overlap for burst
%dlist=44;

%trend
%dlist=115
%dlist=20;

%DATA VG matching MG -56 -105 -20 -64
%dlist=[    44    50        72      84    91   98      112   115 ]
%dlist=[17   45    51         73      86    93   100     113   114  18 65 ]


%for visual latency
%dlist=[  2       4     5             8     9   12      14      20 ...
%44    50    56    64    72    84    91    98   105   112   115]

%noisy LFP
%dlist=[50 56 64 72 84]
%good example trial/trial
%dlist=[14 ]
%dlist=[23 24 28 29]



%vg_rostral
%dlist=[4 6 8 12    14    23    24    28    29    50    64    72 84]
%vg_caudal
%dlist=[ 1 2 3 5 9 20    44    56]
%dlist=[ 1 2  9 20    44    56] %3 5 


%%%%%%%%%%
%all MG
%dlist=[17    18    45    51    57    65    73      86    93   100   107   113   114 ]

%18 65
%less good sessions
%dlist=[57 73 93 107]
%dlist=[18 65]
%dlist=[73      86    93   100   107   113   114 ]

%buildup
%dlist=73;%100
%trend
%dlist=100;


%DATA MG -107 -57
%dlist=[17    18    45    51     65    73      86    93   100     113   114 ]
%DATA MG -107 -57 -18 -65
%dlist=[17   45    51         73      86    93   100     113   114 ]


%RESULTS
%DATA MG -107 -57 +18 +65 for all results
dlist=[17   45    51         73      86    93   100     113   114  18 65 ]


%noisy
%dlist=[18 65 ];

%DATA MG matching VG -57 -107 -18 -65
%dlist=[    45    51           73      86    93   100     113   114 ]


%mg_rostral
%dlist=[25    26    30    31    51    65    73 85 86]
%mg_caudal
%dlist=[17    18    45    57]


%%%%%%%%%%
%data given to Michelle
%dlist=[69 75 88]





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures paper laminar
%examples used
%dlist=[28 64 86]
%dlist=[30]


%figure 2 supp ch 5 and 8 for 64 ch5 for 86
%dlist=[64 86]

%for windows
%033016 5deg ROSTRAL ok  open RF sacc tracking
%dlist=28 %[28 30]%9+[18 19 20 21 22]
%121316 15deg CAUDAL ok
%dlist=9+[47 48 52 54]%[47 48 49 50 51 52 53 54]


%figure CSD
%033016 5deg ROSTRAL ok  open RF sacc tracking
%dlist=9+[19 20]% 21 22]% 21 22];%[18 19 20 21 22]
%all profiles
%dlist=[12      14      20 ...
%   44    50    56    64    72    84    91    98   105   112   115 ...
%17    18    45    51    57    65    73      86    93   100   107   113   114  ]
%23 24 28 29 25 26 30 31 
%increase in SGS layers:113 because of misalignment
%dlist=[112 113]

%figure 1 example
%033016 5deg ROSTRAL ok  open RF sacc tracking
%dlist=28%9+[19 20]%[18 19 20 21 22]
%dlist=9+[21 22]%[18 19 20 21 22]

%figure 6 (build-up/burst)
%dlist=64;
%dlist=100;%MG example should be used



%stats matching sessions
%vg
%dlist=[20  44    50    56    64    72    84    91    98   105   112   115]
%mg
%dlist=[18  45    51    57    65    73    86    93   100   107   113   114 ]


%test latency
%dlist=[14 84 56 44 91 50]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data for Sanjeev V1+Uday's data index
%dlist=[17    18    25    26    30    31    45    51    57    65    73]


%only Uday's data
%dlist=[1:9]


