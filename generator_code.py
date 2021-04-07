import os
import shutil
os.chdir(os.getcwd())


parent_dir=os.getcwd()
os.chdir(parent_dir)
# print(path)
# debug_list = ['CHA.cif']
# list_of_158_no_err_results = ['CHA.cif', 'AVL.cif', 'EWO.cif', 'ITH.cif', 'CSV.cif', 'SBN.cif', 'CAN.cif', 'CZP.cif', 'SZR.cif', 'ITW.cif', 'OSI.cif', 'THO.cif', 'MTT.cif', 'BRE.cif', 'ATS.cif', 'SFF.cif', 'SOD.cif', 'MEP.cif', 'AWW.cif', 'MRT.cif', 'MEL.cif', 'AWO.cif', 'RWR.cif', 'MTN.cif', 'SOS.cif', 'CAS.cif', 'SSY.cif', 'MTF.cif', 'DFT.cif', 'CDO.cif', 'BPH.cif', 'KFI.cif', 'NSI.cif', 'RTH.cif', 'BOF.cif', 'IFR.cif', 'OFF.cif', 'STI.cif', 'CFI.cif', 'MFS.cif', 'AFR.cif', 'ASV.cif', 'AFN.cif', 'TON.cif', 'JRY.cif', 'NAT.cif', 'AEL.cif', 'OBW.cif', 'PHI.cif', 'LTJ.cif', 'UFI.cif', 'SGT.cif', 'EEI.cif', 'SBT.cif', 'RRO.cif', 'ATN.cif', 'MAZ.cif', 'EPI.cif', 'SVV.cif', 'UOZ.cif', 'WEI.cif', 'OKO.cif', 'BOG.cif', 'NAB.cif', 'OSO.cif', 'BEC.cif', 'PON.cif', 'SFH.cif', 'SOF.cif', 'LIO.cif', 'AFS.cif', 'MER.cif', 'ETV.cif', 'BOZ.cif', 'NES.cif', 'GME.cif', 'PUN.cif', 'SFS.cif', 'ATO.cif', 'SAS.cif', 'SAO.cif', 'NPO.cif', 'JSN.cif', 'ZON.cif', 'AFI.cif', 'TER.cif', 'LOV.cif', 'LTA.cif', 'ITE.cif', 'APC.cif', 'SFE.cif', 'GON.cif', 'PWW.cif', 'AFX.cif', 'MVY.cif', 'AEN.cif', 'LTL.cif', 'RWY.cif', 'UEI.cif', 'ATT.cif', 'EDI.cif', 'NPT.cif', 'JSW.cif', 'MEI.cif', 'NON.cif', 'GIS.cif', 'PCR.cif', 'IFY.cif', 'USI.cif', 'LAU.cif', 'ERI.cif', 'BIK.cif', 'SFN.cif', 'ABW.cif', 'APD.cif', 'GOO.cif', 'SAF.cif', 'IHW.cif', '_CON.cif', 'RTE.cif', 'IFO.cif', 'ACO.cif', 'STF.cif', 'EZT.cif', 'JOZ.cif', 'OWE.cif', 'MOR.cif', 'MON.cif', 'BCT.cif', 'AFO.cif', 'JNT.cif', 'ETR.cif', 'DAC.cif', 'AEI.cif', 'UTL.cif', 'CGS.cif', 'YUG.cif', 'PTY.cif', 'SFO.cif', 'AFV.cif', 'AET.cif', 'ANA.cif', 'UOS.cif', 'AHT.cif', 'MTW.cif', 'IWV.cif', 'JBW.cif', 'IWR.cif', 'VSV.cif', 'ATV.cif', 'FER.cif', 'IRR.cif', 'SAV.cif', 'VET.cif', 'DON.cif', 'PWO.cif', 'AST.cif', 'RHO.cif']

for filename in os.listdir(parent_dir):
	if ('.cif' in filename):
		zeolite_name = filename.replace(".cif","")
		path = os.path.join(parent_dir, zeolite_name)
		if os.path.exists(path):
			pass
		else:
			os.mkdir(path)

		file1 = "main_generalized.py"
		file2 = "Opti_Make_matrix_general.py"
		file3 = "Cif_to_Unitcell.py"
		file4 = "convert_xyz_to_csv.py"
		file5 = "indiv_job.sh"
		# file6 = "reader_code.py"
		# file7 = "part2_code.py"
		file8 = "SRU_v2.gms"
		# file9 = "main2.py"

		# os.replace(parent_dir +'/'+ filename, path+'/'+ filename)
		shutil.copyfile(parent_dir +'/'+ filename, path+'/'+ filename)
		shutil.copyfile(parent_dir +'/'+ file1, path+'/'+ file1)
		shutil.copyfile(parent_dir +'/'+ file2, path+'/'+ file2)
		shutil.copyfile(parent_dir +'/'+ file3, path+'/'+ file3)
		shutil.copyfile(parent_dir +'/'+ file4, path+'/'+ file4)
		shutil.copyfile(parent_dir +'/'+ file5, path+'/'+ file5)
		# shutil.copyfile(parent_dir +'/'+ file6, path+'/'+ file6)
		# shutil.copyfile(parent_dir +'/'+ file7, path+'/'+ file7)
		shutil.copyfile(parent_dir +'/'+ file8, path+'/'+ file8)
		# shutil.copyfile(parent_dir +'/'+ file9, path+'/'+ file9)
		
		# os.replace(parent_dir +'/'+ sol_name, path+'/'+ 'hamiltonian_soln.txt')
		# print(path)
		# print(os.getcwd())
		# os.chdir(path)

		# f = open('commandlist','a')
		# f.writelines("cd " + path + "; " + "python main_generalized.py" '\n')
		# f.close()
		# print("cd ", path, "; bsub < indiv_job.sh; cd ..;\n")
		os.system("cd " + str(path) + "; bsub < indiv_job.sh; cd ..;")
