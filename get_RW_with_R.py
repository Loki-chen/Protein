import pandas as pd
import numpy as np
import scipy.sparse as ss
from Protein import get

"""Random walk with restart
   Semiliar neighbor RA
   RA degree RD
"""

def get_matrix_total(data_file):
    # 读取表格
    data = pd.read_excel(data_file)  # 24743 rows  X 2 columns

    # start为表格的第一列
    start = []  # 空列表

    for i in data['start']:  # 此处的data['start'] 是一个名字为start，长度为24743的对象，对象存储的内容为第一列的所有元素。
        if i not in start:
            start.append(i)  # 在列表尾添加对象

    print("the number of start is", len(start))

    # end为表格的第二列
    end = []
    for i in data['end']:
        if i not in end:
            end.append(i)

    print("the number of end is", len(end))

    # total_n为去重后的蛋白质字符串列表
    total_n = []
    # total为原始的蛋白质列表
    total = start + end

    # 列表去重 for循环去重
    for i in total:
        if i not in total_n:
            total_n.append(i)
    # n为total_n的长度

    print("the number of list uniq total is", len(total_n))  # 打印去重后的列表长度

    # 将读取的表格转化为稀疏矩阵
    matrix = ss.lil_matrix((len(total_n), len(total_n)))

    for i in range(len(data['start'])):
        index_start = total_n.index(data['start'][i])
        index_end = total_n.index(data['end'][i])
        matrix[index_start, index_end] = 1  # 给矩阵赋值
        matrix[index_end, index_start] = 1
    return matrix, total_n


def get_Gene12_T(matrix_A, total_n):
    data = pd.DataFrame(pd.read_excel("GENE.xlsx"))
    protein_id = []
    for i in data['protein']:
        protein_id.append(i)
    # print(protein_id)
    count = 0
    for i in range(len(total_n)):
        if total_n[i] in protein_id:
            count += 1
    # print('count', count)
    #  protein ex ex ex....
    list_gene1 = list()
    for i in range(len(protein_id)):
        list_gene2 = list(data.loc[i])
        list_gene1.append(list_gene2)
    # print('data.loc',data.loc[i])
    # print('list_Go1',np.array(list_gene1))
    list_FT = list()
    for i in range(len(protein_id)):
        list_FT1 = list()
        for j in range(1, 13):
            a = list_gene1[i][j]
            b = list_gene1[i][j + 12]
            c = list_gene1[i][j + 24]
            h = (a + b + c) / 3
            list_FT1.append(h)
        list_FT.append(list_FT1)

    AT_list = list()
    for i in range(len(protein_id)):
        u = np.var(list_FT[i][:])
        m = np.mean(list_FT[i][:])
        sigma_noraml = m + 3 * u * (1 - (1 / (1 + pow(u, 2))))
        AT_list.append(sigma_noraml)

    Active_list = list()
    for i in range(len(protein_id)):
        Active_list1 = list()
        for j in range(len(list_FT1)):
            if list_FT[i][j] > AT_list[i]:
                Active_list1.append(list_FT[i][:].index(list_FT[i][j]) + 1)
        Active_list.append(Active_list1)

    # print('Active_list', Active_list)
    Gene_T = dict(zip(protein_id, Active_list))
    # print('Gene_dict', Gene_T.keys())

    index_neighbor = list()
    for i in range(len(total_n)):
        index_list1 = np.argwhere(matrix_A[i,] == 1)
        index_list2 = list()

        for j in index_list1:
            index_list2.append(j[1])
        index_neighbor.append(index_list2)
    # print('index_nei', index_neighbor)

    R_N = ss.lil_matrix((len(total_n), len(total_n)))


    for i in range(len(total_n)):
        for j in index_neighbor[i][:]:
            if total_n[i] in protein_id and total_n[j] in protein_id:
                x = total_n[i]
                y = total_n[j]
                X_index = protein_id.index(x)
                Y_index = protein_id.index(y)
                C = list(set(Active_list[X_index]) & set(Active_list[Y_index]))
                RN = len(C)
                # 属于同一基因周期
                if RN > 0:
                    R_N[i, j] = 1
                    R_N[j, i] = 1
    # print('1-12周期相似下的可靠邻居矩阵', R_N)

    # print(P)

    return R_N


def get_SUBcell_T(matrix_A, total_n):
    loc_size = 11
    locations = ['Cytoskeleton', 'cytosol', 'Endoplasmic', 'endosome', 'extracellular', 'Golgi', 'Mitochondrion',
                 'nucleus',
                 'peroxisome', 'plasma', 'vacuole']
    ix_dict = dict()
    # enumerate:将location列表组合成一个索引序列
    for i, loc in enumerate(locations):
        ix_dict[loc] = i

    # print("ix", ix_dict[loc])
    protein_dict = dict()
    for loc in locations:
        # print('reading %s ...' % loc)
        with open('processed_yeast__(%s).csv' % loc, mode='r') as file:
            file.readline()
            for line in file.readlines():
                protein = line.split(',')[0].strip()
                key = ix_dict[loc]
                if protein in protein_dict:
                    protein_dict[protein][key] = 1
                else:
                    protein_dict[protein] = [0] * loc_size
                    protein_dict[protein][key] = 1

    protein_id = list(protein_dict.keys())

    # print('SUBCell_dict.key', list(protein_dict.keys()))
    SUBcell_T = list(protein_dict.values())
    # print('SUBCell_list', SUBcell_T)

    count = 0
    for i in range(len(total_n)):
        if total_n[i] in protein_dict.keys():
            count += 1
    print(count)

    Sub_T = list()
    for i in range(len(SUBcell_T)):
        Sub1_T = list()
        for j in range(len(SUBcell_T[i])):
            if SUBcell_T[i][j] == 1:
                Sub1_T.append(j)
        Sub_T.append(Sub1_T)

    # print('sub_T', Sub_T)

    index_neighbor = list()
    for i in range(len(total_n)):
        index_list1 = np.argwhere(matrix_A[i,] == 1)
        index_list2 = list()

        for j in index_list1:
            index_list2.append(j[1])
        index_neighbor.append(index_list2)

    R_N1 = ss.lil_matrix((len(total_n), len(total_n)))
    for i in range(len(total_n)):
        for j in index_neighbor[i][:]:
            if total_n[i] in protein_id and total_n[j] in protein_id:
                x = total_n[i]
                y = total_n[j]

                X_index = protein_id.index(x)
                Y_index = protein_id.index(y)

                C = list(set(Sub_T[X_index]) & set(Sub_T[Y_index]))
                RN1 = len(C)
                # 属于同一基因周期
                if RN1 > 0:
                    R_N1[i, j] = 1
                    R_N1[j, i] = 1
    # print('1-12周期相似下的可靠邻居矩阵', R_N1)

    # name = ['Cytoskeleton', 'cytosol', 'Endoplasmic', 'endosome', 'extracellular', 'Golgi', 'Mitochondrion', 'nucleus',
    #         'peroxisome', 'plasma', 'vacuole']
    # pro_id = ['protein']

    # Subfile = pd.DataFrame(columns=name, data=Subcellular_list)
    # Subfile.to_csv("Subcellular.csv")
    #
    # Subfile2 = pd.DataFrame(columns= pro_id, data= Subcellular_idlist)
    # Subfile2.to_csv("Sub_id.csv")
    return R_N1


def get_GO_terms(matrix_A, total_n):
    data = pd.DataFrame(pd.read_excel("Go_data.xlsx"))
    Go_dict = dict()
    protein_id = []
    for i in data['protein']:
        protein_id.append(i)
    # print('protein_IDnex',protein_id)
    count = 0
    for i in range(len(total_n)):
        if total_n[i] in protein_id:
            count += 1
    # print('count', count)
    #  protein ex ex ex....
    list_Go1 = list()
    for i in range(len(protein_id)):
        list_Go2 = list(data.loc[i][1:35])
        list_Go1.append(list_Go2)
    # print('data.loc',len(data.loc[i]))
    # print('list_Go1',np.array(list_Go1))

    A = np.array(list_Go1)
    B = np.nan_to_num(A)

    Go_term = list()
    for i in range(len(B)):
        Go_term1 = list()
        for j in range(len(B[i])):
            if B[i][j] != 0:
                Go_term1.append(B[i][j])
        Go_term.append(Go_term1)
    # print('list_GO_remove nan', Go_term)

    index_neighbor = list()
    for i in range(len(total_n)):
        index_list1 = np.argwhere(matrix_A[i,] == 1)
        index_list2 = list()

        for j in index_list1:
            index_list2.append(j[1])

        index_neighbor.append(index_list2)
    # print('index_nei', index_neighbor)
    GO_list = list()
    for i in range(len(total_n)):
        GO_list1 = list()
        for j in index_neighbor[i][:]:
            if total_n[i] in protein_id and total_n[j] in protein_id:
                x = total_n[i]
                y = total_n[j]
                X_index = protein_id.index(x)
                Y_index = protein_id.index(y)
                A, B = len(Go_term[X_index]), len(Go_term[Y_index])
                if A > 0 and B > 0:
                    C = len(list(set(Go_term[X_index]) & set(Go_term[Y_index])))
                    P = (C ** 2)  / (A * B)
                    GO_list1.append(P)
                else:
                    GO_list1.append(0)
            else:
                P = 0
                GO_list1.append(P)
        GO_list.append(GO_list1)

    GO = list()
    for i in range(len(total_n)):
        b = GO_list[i]
        GO.append(sum(b))
    print('Go_list', GO)
    return GO


def get_SCS(total_n):

    loc_size = 11
    locations = ['Cytoskeleton', 'cytosol', 'Endoplasmic', 'endosome', 'extracellular', 'Golgi', 'Mitochondrion',
                 'nucleus',
                 'peroxisome', 'plasma', 'vacuole']
    ix_dict = dict()
    # enumerate:将location列表组合成一个索引序列
    for i, loc in enumerate(locations):
        ix_dict[loc] = i
    # print("ix", ix_dict[loc])
    protein_dict = dict()
    for loc in locations:
        # print('reading %s ...' % loc)
        with open('processed_yeast__(%s).csv' % loc, mode='r') as file:
            file.readline()
            for line in file.readlines():
                protein = line.split(',')[0].strip()
                key = ix_dict[loc]
                if protein in protein_dict:
                    protein_dict[protein][key] = 1
                else:
                    protein_dict[protein] = [0] * loc_size
                    protein_dict[protein][key] = 1

    Subcellular_idlist = list(protein_dict.keys())
    Subcellular_list = list(protein_dict.values())
    # print("Subcellular_list:", Subcellular_list)

    # print('protein_dict', protein_dict)

    # name = ['Cytoskeleton', 'cytosol', 'Endoplasmic', 'endosome', 'extracellular', 'Golgi', 'Mitochondrion', 'nucleus',
    #        'peroxisome', 'plasma', 'vacuole']
    # Subfile = pd.DataFrame(columns=name, data=Subcellular_list)
    # Subfile.to_csv("Subcellular.csv")
    SN_list1 = pd.read_excel('Cytoskeleton.xlsx')
    SN_list2 = pd.read_excel('cytosol.xlsx')
    SN_list3 = pd.read_excel('Endoplasimic.xlsx')
    SN_list4 = pd.read_excel('endosome.xlsx')
    SN_list5 = pd.read_excel('extracellular.xlsx')
    SN_list6 = pd.read_excel('Golgi.xlsx')
    SN_list7 = pd.read_excel('Mitochondrion.xlsx')
    SN_list8 = pd.read_excel('nucleus.xlsx')
    SN_list9 = pd.read_excel('peroxisome.xlsx')
    SN_list10 = pd.read_excel('plasma.xlsx')
    SN_list11 = pd.read_excel('vacuole.xlsx')

    Cyt_list = list(SN_list1['Q0045'])
    cyt_list = list(SN_list2['Q0045'])
    End_list = list(SN_list3['Q0045'])
    end_list = list(SN_list4['Q0045'])
    ext_list = list(SN_list5['Q0045'])
    Gol_list = list(SN_list6['Q0045'])
    Mit_list = list(SN_list7['Q0045'])
    nuc_list = list(SN_list8['Q0045'])
    per_list = list(SN_list9['Q0045'])
    pla_list = list(SN_list10['Q0045'])
    vac_list = list(SN_list11['Q0045'])

    Cytoskeleton_list = list()
    cytosol_list = list()
    Endoplasmic_list = list()
    endosome_list = list()
    extracellular_list = list()
    Golgi_list = list()
    Mitochondrion_list = list()
    nucleus_list = list()
    peroxisome_list = list()
    plasma_list = list()
    vacuole_list = list()
    SN_max_list = list()
    # 每个亚细胞中出现的蛋白质数量
    for i in total_n:
        if i in Cyt_list:
            Cytoskeleton_list.append(i)
        if i in cyt_list:
            cytosol_list.append(i)
        if i in End_list:
            Endoplasmic_list.append(i)
        if i in end_list:
            endosome_list.append(i)
        if i in ext_list:
            extracellular_list.append(i)
        if i in Gol_list:
            Golgi_list.append(i)
        if i in Mit_list:
            Mitochondrion_list.append(i)
        if i in nuc_list:
            nucleus_list.append(i)
        if i in per_list:
            peroxisome_list.append(i)
        if i in pla_list:
            plasma_list.append(i)
        if i in vac_list:
            vacuole_list.append(i)
        if i in Subcellular_idlist:
            SN_max_list.append(i)

    # # 子细胞定位指标：

    Score_Cyt = list()
    Score_cyt = list()
    Score_End = list()
    Score_end = list()
    Score_ext = list()
    Score_Gol = list()
    Score_Mit = list()
    Score_nuc = list()
    Score_per = list()
    Score_pla = list()
    Score_vac = list()

    for i in total_n:
        if i in Cyt_list:
            Score_Cyt.append(len(Cytoskeleton_list) / len(SN_max_list))
        else:
            Score_Cyt.append(0)

        if i in cyt_list:
            Score_cyt.append(len(cytosol_list) / len(SN_max_list))
        else:
            Score_cyt.append(0)

        if i in End_list:
            Score_End.append(len(Endoplasmic_list) / len(SN_max_list))
        else:
            Score_End.append(0)

        if i in end_list:
            Score_end.append(len(endosome_list) / len(SN_max_list))
        else:
            Score_end.append(0)

        if i in ext_list:
            Score_ext.append(len(extracellular_list) / len(SN_max_list))
        else:
            Score_ext.append(0)

        if i in Gol_list:
            Score_Gol.append(len(Golgi_list) / len(SN_max_list))
        else:
            Score_Gol.append(0)

        if i in Mit_list:
            Score_Mit.append(len(Mitochondrion_list) / len(SN_max_list))
        else:
            Score_Mit.append(0)

        if i in nuc_list:
            Score_nuc.append(len(nucleus_list) / len(SN_max_list))
        else:
            Score_nuc.append(0)

        if i in per_list:
            Score_per.append(len(peroxisome_list) / len(SN_max_list))
        else:
            Score_per.append(0)

        if i in pla_list:
            Score_pla.append(len(plasma_list) / len(SN_max_list))
        else:
            Score_pla.append(0)

        if i in vac_list:
            Score_vac.append(len(vacuole_list) / len(SN_max_list))
        else:
            Score_vac.append(0)

    Sum_list = [a + b + c + d + e + f + g + h + i + j + k for a, b, c, d, e, f, g, h, i, j, k in zip(
        Score_Cyt, Score_cyt, Score_End, Score_end, Score_ext, Score_Gol, Score_Mit,
        Score_nuc, Score_per, Score_pla, Score_vac)]
    print('SCS', Sum_list)
    return Sum_list


# 计算转移概率矩阵
def Probability_matrix(matrix_A, total_n):

    # 相似周期邻接矩阵
    R_AN = ss.lil_matrix((len(total_n), len(total_n)))
    R_N = get_Gene12_T(matrix_A, total_n)
    R_N1 = get_SUBcell_T(matrix_A, total_n)

    for i in range(len(total_n)):
        for j in range(len(total_n)):
            if R_N[i, j] == 1 and R_N1[i, j] == 1:
                R_AN[i, j] = 3 * matrix_A[i, j]
                R_AN[j, i] = 3 * matrix_A[j, i]
            elif R_N[i, j] == 1 or R_N1[i, j] == 1:
                R_AN[i, j] = 2 * matrix_A[i, j]
                R_AN[j, i] = 2 * + matrix_A[j, i]
            else:
                R_AN[i, j] = matrix_A[i, j]
                R_AN[j, i] = matrix_A[j, i]
        R_AN[i, i] = 1
    # print('相似邻居矩阵', R_AN)



    R_spend = list()
    for i in range(R_AN.shape[0]):
        R_spend.append(np.sum(R_AN[i, ]))

    R_D = ss.lil_matrix((len(total_n), len(total_n)))
    for i in range(len(total_n)):
            R_D[i, i] = 1 / (R_spend[i] + 1)



    # print('相似度矩阵', R_D)
    # 转移概率矩阵
    P = np.dot(R_AN, R_D).todense()

    # print('P', np.array(P))
    print('Matrix-done!')
    return P


def get_RW_with_R(a, P_matrix, matrix_A, total_n):

    SCS = get_SCS(total_n)
    GOS = get_GO_terms(matrix_A, total_n)

    P = list()
    count = 0
    for i in range(len(total_n)):
        P1 = list()
        for j in range(len(total_n)):
            P1.append(P_matrix[i, j])
            if P_matrix[i, j] != 0:
                count += 1
        P.append(P1)

    INI = list()
    for i in range(len(total_n)):
        INI.append(SCS[i] * GOS[i])

    RWR = [1 for index in range(len(total_n))]

    b = 0.000001
    Times = 0
    Conditional = np.inf
    while Conditional > b:
        A = np.dot(P, RWR)  # N * 1 dim
        R_T = RWR
        RWR = (1 - a) * np.array(A) + a * np.array(INI)
        R_T1 = RWR
        D = np.array(R_T1) - np.array(R_T)
        Conditional = np.linalg.norm(D)
        Times += 1
    print('迭代次数', Times)

    return RWR


def get_RWR_reverse(t):

    total_index = dict()
    for i in range(len(t)):
        total_index[i] = t[i]
    SNDC_reverse_list = sorted(total_index.items(), key=lambda asd: asd[1], reverse=True)
    return SNDC_reverse_list


def Test(reverse_list, matrix, total_n):
    count = int(len(total_n) * 0.01)
    # count = 100
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")
    count = int(len(total_n) * 0.05)
    # count = 200
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")
    count = int(len(total_n) * 0.1)
    # count = 300
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")
    count = int(len(total_n) * 0.15)
    # count = 400
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")
    count = int(len(total_n) * 0.20)
    # count = 500
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")
    count = int(len(total_n) * 0.25)
    # count = 600
    ess_count, no_ess_count, ess_count1, no_ess_count1 = get.get_ess_count(reverse_list, total_n, count)
    print(count, "个蛋白质有", ess_count, "个必须蛋白质,有", no_ess_count, "非必须蛋白质")





if __name__ == "__main__":
    # data_1 = "Krogan.xlsx"   # 363  0.6   1110.0      SAS
#     data_2 = "YDIP.xlsx"
#     data_3 = "YMIPS.xlsx"   # 451  1   1110.0
#     data_4 = 'Gavin.xlsx'   # 258  0.6  1110.0         366  600  1     !!!!!!!!!        SAS
# #     data_5 = 'BioGRID.xlsx'   # 353.  1 600 !!!!!!!!!!          SAS
      data_6 = 'YHQ.xlsx'
#
      x, y = get_matrix_total(data_6)
      p = Probability_matrix(x, y)
#     # for a in range(0, 11, 1):
      t = get_RW_with_R(0.5, p, x, y)
      Reverse = get_RWR_reverse(t)
      print("%s网络第%f次的输出结果########################################"%(data_6, 0.5))
      Test(Reverse, x, y)
    #     print("%s网络第%f次的输出结果########################################"%(data_5, a/10))



