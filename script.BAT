
 
set list=10	11	12	13	14	15	16	17	18	19	1	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	3	40	41	42	43	44	45	4	54	55	56	5	60	6	7	8	9
set list=28
(for %%a in (%list%) do ( 

	C:\\Users\\liors\\source\\repos\\ctrm\\x64\\Release\\ctrm.exe 60 1 30 -40 20 39 1 1 -1 -1 100 20 100 4 4 C:\\Users\\liors\\source\\repos\\ctrm\\ctrm\\aquisition38chan.txt C:\\Users\\liors\\source\\repos\\ctrm\\ctrm\\vrsMemadion.txt C:\\Users\\liors\\source\\repos\\ctrm\\ctrm\\MemadionData\\LB_28N ent

	move C:\\Users\\liors\\source\\repos\\ctrm\\ctrm\\MemadionData\\*.npy results
))
