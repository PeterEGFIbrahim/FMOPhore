###################################################################################################
# Analysis of the Complex 
###################################################################################################
class ComplexAnalysis:
    def __init__(self, pdb_file, same_target=False):
        self.pdb_file = pdb_file
        self.same_target = same_target

    def com_analyze(self):
        import numpy as np
        import seaborn as sns
        import pandas as pd
        from itertools import islice
        from math import factorial
        import matplotlib.pyplot as plt
        import re
        import shutil
        import csv
        from numpy import nan
        import glob
        import os

        output_lines = [] 
        inp_file_path = [f'{self.pdb_file[:-4]}.inp']
        if inp_file_path:
            inp_file_path = inp_file_path[0]
        else:
            print("No .inp file found in the current directory.")
            exit(1)
        with open(inp_file_path, 'r') as file:
            nfrag_found = False
            for line in file:
                if re.search(r'\bNFRAG\b', line):
                    nfrag_found = True
                    output_lines.append(line.rstrip()) 
                    break  
            if nfrag_found:
                inside_block = True
                for line in file:
                    if re.search(r'\bINDAT\b', line):
                        output_lines.append(line.rstrip())
                        break
                    if inside_block and line.strip():
                        line = re.sub(r'\bLIG-\d+\b', 'LIG', line)
                        output_lines.append(line.rstrip())
                    if re.search(r'\bNFRAG\b', line):
                        inside_block = True
            else:
                output_lines.append("NFRAG not found in the file.")
        log_file_path = glob.glob(f'{self.pdb_file[:-4]}.log')
        if log_file_path:
            log_file_path = log_file_path[0]
        else:
            print("No .log file found in the current directory.")
            exit(1)
        with open(log_file_path, 'r') as existing_file:
            existing_lines = existing_file.readlines()
            for i in range(len(existing_lines)):
                existing_lines[i] = re.sub(r'\bLIG-\d+\b', 'LIG', existing_lines[i])
        with open(log_file_path, 'w') as output_file:
            for line in output_lines:
                output_file.write(line + '\n')
            output_file.writelines(existing_lines)
        output_file='ligand_heatmap.csv'
        FMO_file=[f'{self.pdb_file[:-4]}.log']
        selection='LIG'
        total_list=[]
        Ees_list=[]
        Eex_list=[]
        Ect_list=[]
        Edisp_list=[]
        Gsol_list=[]
        for file in FMO_file:
            pattern_01 = re.compile(r'(\s{4}I\s{4}J\sDL)')
            pattern_02 = re.compile(r'(\s{6}NFRAG=\d+)')     
            with open(output_file,'w+') as f:
                for i, line in enumerate(open(file)): 
                    for match in re.finditer(pattern_01, line):
                        start_line=(i+3)
                    for match in re.finditer(pattern_02, line):
                        bla=match.groups()
                        nfrag=bla[0].split("=")[1]            
                        nlines=(factorial(int(nfrag))/(factorial(2)*factorial(int(nfrag)-2)  ))      
            f.close
            with open(output_file,'w+') as f:
                for i, line in enumerate(open(file)): 
                    if i in range(start_line-1, int(start_line+nlines)):
                                f.write(line)
            f.close
            pattern_01 = re.compile(r'(\s{6}FRGNAM\(\d\)=)')
            pattern_02 = re.compile(r'(\s{6}INDAT)')  
            pattern_03 = re.compile(r'(\s{6})')  
            combined_pat = [r'(\s{6}FRGNAM\(\d\)=)',
                            r'(\s{6}INDAT)',
                            r'(\s{6})']
            residues=[]
            for i, line in enumerate(open(file)):
                for match in re.finditer(pattern_01, line):
                    start_line=i
                for match in re.finditer(pattern_02, line):
                    end_line=i
            f.close
            for i, line in enumerate(open(file)):

                if i in range(start_line, end_line):
                    residues.append(re.sub('|'.join(combined_pat)," ", line))
            f.close 
        ###############################################################################
            residues=(list(map(lambda x:x.strip(),residues)))
            residues = [x.strip(' ') for x in residues]
            my_residues=[]
            for element in residues:
                my_residues.extend(element.split(','))
            my_residues = [x.strip(' ') for x in my_residues]
            my_residues = [x for x in my_residues if x != '']
            fragment_num=str(my_residues.index(selection)+1)
            df = pd.DataFrame()
            for i, line in enumerate(open(output_file)):
                if fragment_num in (line[:10]):
                    temp=pd.DataFrame(
                        {
                            "I":[line[:5]],
                            "J":[line[6:10]],
                            "DL":[line[11:13]],
                            "Z":[line[14:16]],
                            "R":[line[17:23]],
                            "Q(I->J)":[line[24:31]],
                            "EIJ-EI-EJ":[line[32:41]],
                            "dDIJ*VIJ":[line[42:50]],
                            "total":[line[51:60]],
                            "Ees":[line[61:70]],
                            "Eex":[line[71:79]],
                            "Ect+mix":[line[80:88]],
                            "Edisp":[line[89:97]],
                            "Gsol":[line[98:106]]

                        }
                    )
                    df = pd.concat([df, temp])
            my_residues.pop(int(fragment_num)-1)
            df.insert(0,"Fragment",my_residues)
            df_map=df[['Fragment','total']].T
            df_map.columns = df_map.iloc[0]
            df_map=df_map.drop(['Fragment'],axis=0)
            df_map = df_map=df_map.set_index([pd.Index([(os.path.splitext(file)[0])])])
            df_map = df_map.apply(pd.to_numeric, errors='coerce')
            if 'ligands' not in locals():
                ligands=[df_map.index[0]]
            else:
                ligands.append(df_map.index[0])
            if len(ligands)==1:
                df_map_total=df_map
            else:
                df_map_total=df_map_total.append(df_map,sort='false') #may or may not need sort here!
                df_map_total=df_map_total.reindex(sorted(df_map_total.columns, key=lambda s: s[s.rfind('-'):]), axis=1)
                df_map_total=df_map_total.reindex(sorted(df_map_total.columns, key=lambda x: int(re.search(r'\d+$',x).group())), axis=1)
                df_map_total.index=df_map_total.index.str.replace('_FMO','')
                df_map_total.index=df_map_total.index.str.replace('min','')
            cols=df.columns.drop(['Fragment','DL'])
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
            df2=df
            sorter=list(df_map_total.columns) 
            true_sort = [s for s in sorter if s in df2.Fragment.unique()]
            df2=df2.set_index('Fragment').loc[true_sort].reset_index()
            df2.set_index("Fragment",drop=True,inplace=True)  
            if 'Ees' and 'Eex' and 'Ect+mix' and 'Edisp' and 'Gsol' not in locals():
                Ees=[list(df2['Ees'])]
                Eex=[list(df2['Eex'])]
                Ect=[list(df2['Ect+mix'])]
                Edisp=[list(df2['Edisp'])]
                Gsol=[list(df2['Gsol'])]
            else:
                Ees.append(list(df2['Ees']))
                Eex.append(list(df2['Eex']))
                Ect.append(list(df2['Ect+mix']))
                Edisp.append(list(df2['Edisp'])) 
                Gsol.append(list(df2['Gsol']))
            total_list.append(df['total'].sum())
            Ees_list.append(df['Ees'].sum())
            Eex_list.append(df['Eex'].sum())
            Ect_list.append(df['Ect+mix'].sum())
            Edisp_list.append(df['Edisp'].sum())
            Gsol_list.append(df['Gsol'].sum())
        bla=list(map(tuple, np.where(np.isnan(df_map_total.values.tolist()))))
        bla_a=bla[0]
        bla_b=bla[1]
        for f, b in zip(bla_a, bla_b):
            Ees[f].insert(b, nan)  
            Eex[f].insert(b, nan) 
            Ect[f].insert(b, nan) 
            Edisp[f].insert(b, nan) 
        for file in FMO_file:
            pattern_01 = re.compile(r'(\s{4}I\s{4}J\sDL)')
            pattern_02 = re.compile(r'(\s{6}NFRAG=\d+)')  
            with open(output_file,'w+') as f:
                for i, line in enumerate(open(file)):
                    for match in re.finditer(pattern_01, line):
                        start_line=(i+3)

                    for match in re.finditer(pattern_02, line):
                        bla=match.groups()
                        nfrag=bla[0].split("=")[1]
                        nlines=(factorial(int(nfrag))/(factorial(2)*factorial(int(nfrag)-2) ))
            f.close
            with open(output_file,'w+') as f:
                for i, line in enumerate(open(file)):
                    if i in range(start_line-1, int(start_line+nlines)):
                        f.write(line)
            f.close
            pattern_01 = re.compile(r'(\sINPUT CARD\>\s{6}FRGNAM\(\d\)=)')
            pattern_02 = re.compile(r'(\sINPUT CARD\>\s{6}INDAT)')  
            pattern_03 = re.compile(r'(\sINPUT CARD\>\s{6})')
            combined_pat = [
                r'(\sINPUT CARD\>\s{6}FRGNAM\(\d\)=)',
                r'(\sINPUT CARD\>\s{6}INDAT)',
                r'(\sINPUT CARD\>\s{6})'
            ]
            residues=[]
            for i, line in enumerate(open(file)):
                for match in re.finditer(pattern_01, line):
                    start_line=i
                for match in re.finditer(pattern_02, line):
                    end_line=i
            f.close
            for i, line in enumerate(open(file)):
                if i in range(start_line, end_line):
                    residues.append(re.sub('|'.join(combined_pat)," ", line))
            f.close
        df_map_total.insert(0,'total', total_list)
        df_map_total=df_map_total.sort_values(by=['total'])
        df_map_total
        df_map_total.to_csv (r'export_dataframe.csv', header=True)
        df_map_total
        selected = df_map_total[df_map_total["total"] >= -400]
        selected=selected.drop(columns=['total'])
        selected
        cmap_colors = ['#FF0000', '#FF8800', '#FFFF00', '#00FF00', '#0000FF']
        plt.style.use('seaborn-darkgrid')
        plt.figure(figsize=(20,10))
        ax = sns.heatmap(selected, linewidths=.5, cmap="RdYlGn", vmin=-25, vmax=5, center=0)
        cbar = ax.collections[0].colorbar
        cbar.set_label("kcal/mol")
        plt.tick_params(labelsize=12)
        plt.yticks(rotation='horizontal')
        plt.xlabel('')
        plt.tight_layout(pad=2)
        plt.savefig('DFTB_heatmap_lig.png')
        plt.show()
        ###################################################################################################
        cols=df.columns.drop(['Fragment','DL'])
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
        df2=(df.loc[df['total'] <= 50])
        df2.set_index("Fragment",drop=True,inplace=True)

        ###################################################################################################
        import matplotlib.pyplot as plt
        import os
        # df2 = df2.sort_values(by=['total'])
        df2.iloc[:, 9:14].plot.bar(stacked=True, color=('gold', 'mediumseagreen', 'firebrick', 'steelblue', 'darksalmon'))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1)
        plt.xlabel("Binding residues")
        plt.ylabel("kcal/mol")
        plt.tight_layout(pad=2)
        plt.savefig((FMO_file)[0] + "_PIEDA.png")
        df3 = df2
        df3.iloc[:, 8:9].plot.bar(color=('darkred'))
        plt.xlabel("Binding residues")
        plt.ylabel("kcal/mol")
        plt.tight_layout(pad=2)
        plt.savefig((os.path.splitext(FMO_file[0])[0]) + "_total.png")
        plt.show()
        df.to_csv('PIEDA_values.txt', sep='\t',float_format='%.3f')
        # ###################################################################################################
        # Transpose the data
        # ###################################################################################################
        df_subset = df.iloc[:, [0] + list(range(9, 15))].reset_index(drop=True)
        df_subset.to_csv('Binding_site_energies.txt', sep='\t', float_format='%.3f', index=False)
        ###################################################################################################
        ###################################################################################################
        fragments=list((range(1,len(my_residues)+2)))
        fragments.remove(int(fragment_num))
        total_list.append(df['total'].sum())
        Ees_list.append(df['Ees'].sum())
        Eex_list.append(df['Eex'].sum())
        Ect_list.append(df['Ect+mix'].sum())
        Edisp_list.append(df['Edisp'].sum())
        Gsol_list.append(df['Gsol'].sum())
        with open("totals.csv","w") as f:
            f.write("index,DFTB_total,Ees_DFTB,Eex_DFTB,Ect_DFTB,Edisp_DFTB,Gsol_DFTB\n")
            for (index,DFTB_total,Ees_DFTB,Eex_DFTB,Ect_DFTB,Edisp_DFTB,Gsol_DFTB) in zip(df_map_total.index,total_list,Ees_list,Eex_list,Ect_list,Edisp_list,Gsol_list):
                f.write("{0},{1:.2f},{2:.2f},{3:.2f},{4:.2f},{5:.2f},{6:.2f}\n".format(index,DFTB_total,Ees_DFTB,Eex_DFTB,Ect_DFTB,Edisp_DFTB,Gsol_DFTB))
        df2 = df_map_total[["total"]].mean()
        with open("total_average.txt", "a") as outfile:
            outfile.write("Total_average {} \n".format(str(df2)))
        # ###################################################################################################
        # # Analysis of binding energy
        # ###################################################################################################
        with open(log_file_path, 'r') as existing_file:
            existing_lines = existing_file.readlines()
        one_body = None
        two_body = None
        one_body_gsol = None
        two_body_gsol = None
        two_body_interaction = None
        Total_EunDsol2 = None
        for line in existing_lines:
            if 'Total energy of the molecule: Eunco+D(1)=' in line:
                one_body = float(line.split()[-1])
            elif 'Total energy of the molecule: EunD+so(2)=' in line:
                Total_EunDsol2 = float(line.split()[-1])
            elif 'Total energy of the molecule: Eunco+D(2)=' in line:
                two_body = float(line.split()[-1])
            elif 'Total Gsol(1)=' in line:
                one_body_gsol = float(line.split()[-2])
            elif 'Total Gsol(2)=' in line:
                two_body_gsol = float(line.split()[-2])
            elif 'Total interaction (PL state)         E\'int' in line:
                two_body_interaction = float(line.split()[-1])
        if Total_EunDsol2 is not None and one_body is not None and two_body is not None and one_body_gsol is not None and two_body_gsol is not None and two_body_interaction is not None:
            complex_energy = (one_body * 627.51) + one_body_gsol + two_body_interaction
            with open('complex_energy.txt', 'w') as f:
                f.write(f'complex_energy: {complex_energy} kcal/mol')
        else:
            print('Error: Unable to calculate total energy.')
        ###################################################################################################
        # Analysis of Interactions for library plot
        ###################################################################################################
        if self.same_target:
            new_directory = "../../../library_analysis"
            os.makedirs(new_directory, exist_ok=True)
            new_log_file_path = os.path.join(new_directory, os.path.basename(log_file_path))
            new_output_file_path = os.path.join(new_directory, os.path.basename(output_file))
            shutil.copy(log_file_path, new_log_file_path)
            shutil.copy(output_file, new_output_file_path)
        else:
            pass

if __name__ == "__main__":
    """
    FMOPhore V 0.1 - ComplexAnalysis - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = ComplexAnalysis(self.pdb_file)
    pdb_processor.com_analyze()
