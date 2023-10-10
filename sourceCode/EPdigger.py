import numpy as np
import os
from Bio import SeqIO
import math
import re
from PySide2.QtWidgets import QApplication, QMessageBox, QFileDialog
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import QFile, Signal, QObject
from threading import Thread
import os.path
import openpyxl
from PySide2.QtGui import QIcon


# 信号库
class SignalStore(QObject):
    # 定义一种信号
    progress_update = Signal(int)
    # 线程完成信号
    finished = Signal(list)
    # 出错信号类别
    error_window = Signal(int)


# 还可以定义其他作用的信号


# 实例化信号库
so = SignalStore()


class AboutWidget():
    def __init__(self):
        # 从文件加载UI定义
        ui_file = QFile('aboutWin.ui')  # 调用UI文件
        ui_file.open(QFile.ReadOnly)
        ui_file.close()

        # 从 UI 定义中动态 创建一个相应的窗口对象
        # 注意：里面的控件对象也成为窗口对象的属性了
        # 比如 self.ui.button , self.ui.textEdit
        self.ui = QUiLoader().load(ui_file)


class MainWindow():

    def __init__(self):
        # 连接信号到处理的slot函数
        so.progress_update.connect(self.setProgress)
        so.finished.connect(self.show_message)
        so.error_window.connect(self.window_error)

        # 从文件加载UI定义
        ui_file = QFile('mainWin.ui')  # 调用UI文件
        ui_file.open(QFile.ReadOnly)
        ui_file.close()

        # 从 UI 定义中动态 创建一个相应的窗口对象
        # 注意：里面的控件对象也成为窗口对象的属性了
        # 比如 self.ui.button , self.ui.textEdit
        self.ui = QUiLoader().load(ui_file)

        # 文件路径存放
        self.filePaths_ref = []
        self.filePaths_alfa = []
        self.filePaths_rsa = []

        #self.ui.lineEdit.setReadOnly(True)
        self.ui.choice_Button_fa.clicked.connect(self.get_path_ref)

        #self.ui.lineEdit_2.setReadOnly(True)
        self.ui.choice_Button_alfa.clicked.connect(self.get_path_alfa)

        #self.ui.lineEdit_3.setReadOnly(True)
        self.ui.choice_Button_rsa.clicked.connect(self.get_path_rsa)

        # 进度是 0 - 10，
        self.ui.progressBar.setRange(0, 10)

        # run按钮
        self.ui.runButton.clicked.connect(self.runp)

        # 跳转到About窗口
        self.ui.actionAbout.triggered.connect(self.open_AboutWidget)

    def open_AboutWidget(self):
        # 实例化另外一个窗口
        self.aboutw = AboutWidget()
        # 显示新窗口
        self.aboutw.ui.show()

    # 处理进度的slot函数
    def setProgress(self, value):
        self.ui.progressBar.setValue(value)

    # 处理出错提示的slot函数
    def window_error(self, number):
        # 弹窗提示全部等长完整序列出错，请仔细检查
        message = ""
        ew = QMessageBox()
        ew.setWindowTitle("Error")
        ew.setIcon(QMessageBox.Warning)
        if number == 1:
            message = "Error in all full-length sequence file"
        elif number == 2:
            message = "Error in RSA file"
        elif number == 3:
            message = "Please check that these files are correct"
        elif number == 4:
            message = "All full-length sequence file have sequences containing unknown or non-standard amino acids"
        elif number == 5:
            message = "Error in All full-length sequences have only one sequence"
        elif number == 6:
            message = "Please close pv_ASA.xlsx"
        elif number == 7:
            message = "Please close w_value.xlsx"
        elif number == 8:
            message = "Please close result.xlsx"
        else:
            message = "Error in reference sequence file"
        ew.setText(message)
        ew.exec_()

    # 程序运行完毕结果窗口
    def show_message(self, message):
        # 创建一个消息框对象，并设置标题，文本和图标
        result_msg = ""
        msg = QMessageBox()
        msg.setWindowTitle("Result")
        for i in range(len(message[0])):
            result_msg += (message[0][i] + '\n' + message[1][i]) + '\n'
        msg.setText(result_msg)
        msg.setIcon(QMessageBox.Information)
        # 显示消息框
        msg.exec_()

    def get_path_ref(self):
        path1, _ = QFileDialog.getOpenFileName(self.ui, "Select file", "", "Fasta Files (*.Fasta)")
        self.filePaths_ref.append(path1)
        filename = os.path.basename(path1)
        self.ui.label_5.setText(filename)


    def get_path_alfa(self):
        path2, _ = QFileDialog.getOpenFileName(self.ui, "Select file", "", "Fasta Files (*.Fasta)")
        self.filePaths_alfa.append(path2)
        filename = os.path.basename(path2)
        self.ui.label_6.setText(filename)


    def get_path_rsa(self):
        path3, _ = QFileDialog.getOpenFileName(self.ui, "Select file", "", "Text Files (*.txt)")
        self.filePaths_rsa.append(path3)
        filename = os.path.basename(path3)
        self.ui.label_7.setText(filename)


    def runp(self):
        # 窗口长度
        ckLengt = self.ui.spinBox.value()

        # 窗口信号数值

        # 定义二十维坐标
        def d20(seq):
            # 字典存储氨基酸对应的矩阵
            aa_dict = {
                'A': [1.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'C': [0, 1.002, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'D': [0, 0, 1.003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'E': [0, 0, 0, 1.004, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'F': [0, 0, 0, 0, 1.005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'G': [0, 0, 0, 0, 0, 1.006, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'H': [0, 0, 0, 0, 0, 0, 1.007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'I': [0, 0, 0, 0, 0, 0, 0, 1.008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'K': [0, 0, 0, 0, 0, 0, 0, 0, 1.009, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'L': [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.010, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'M': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.011, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'N': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.012, 0, 0, 0, 0, 0, 0, 0, 0],
                'P': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.013, 0, 0, 0, 0, 0, 0, 0],
                'Q': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.014, 0, 0, 0, 0, 0, 0],
                'R': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.015, 0, 0, 0, 0, 0],
                'S': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.016, 0, 0, 0, 0],
                'T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.017, 0, 0, 0],
                'V': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.018, 0, 0],
                'W': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.019, 0],
                'Y': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.020]}

            # 将序列中的每个氨基酸与字典相关联
            result = [aa_dict[aa] for aa in seq]

            # 将列表逐个相加最终得到20D坐标
            a = 0
            b = 0
            lst1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            lst2 = []
            while a < len(result):
                lst1 = list(np.add(lst1, result[b]))
                a += 1
                b += 1

                lst2.append(lst1)  # 得到了各个点的坐标位置

            # 把各个点的坐标相加，方便用于计算pv
            a1 = 0
            b1 = 0
            lst3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            while a1 < len(lst2):
                lst3 = list(np.add(lst3, lst2[b1]))
                a1 += 1
                b1 += 1

            return lst3

        # 图形滑动窗口，计算PV
        def gswnPv(seq):
            lt = len(seq)
            n = ckLengt  # 蛋白质链中氨基酸的数量
            prl = []  # 存放pr
            for a in range(lt):  # 不能大于序列的长度,这里a为下标从0开始
                # 取窗口大小为调节器所得的数
                left = a  # 索引下标左
                right = a + (ckLengt - 1)  # 索引下标右
                if right > lt - 1:
                    # print("窗口扫描结束")
                    break
                elif right < (ckLengt - 1):
                    # print("序列长度小于12")
                    break
                seq1 = seq[left:right + 1]
                seq2 = d20(seq1)  # 返回值为一个抽象曲线上各个点的坐标相加的和
                ui = [(x / n) ** 2 for x in seq2]
                pr = math.sqrt(sum(ui))
                prl.append(pr)
            return prl

        # 1/PV
        def re1pv(pv):
            pv1 = []
            for i in pv:
                p = 1 / i
                pv1.append(p)
            return pv1

        # 归一化公式 y = (x - x_min) / (x_max - x_min) * 10
        def nor(pv1):
            c = max(pv1)
            d = min(pv1)
            b = c - d
            ynl = []
            if b != 0:
                for x in pv1:
                    a = x - d
                    yn = a / b * 10
                    ynl.append(yn)
            else:
                ynl = [1 for i in pv1]

            return ynl

        def gswnAsa(seq):
            lt = len(seq)
            asal = []  # 存放asa

            for a in range(lt):  # 不能大于序列的长度
                # 取窗口大小为12
                left = a  # 索引下标左
                right = a + (ckLengt - 1)  # 索引下标右
                if right > lt - 1:
                    # print("窗口扫描结束")
                    break
                elif right < (ckLengt - 1):
                    # print("序列长度小于12")
                    break
                seq1 = seq[left:right + 1]
                x = 0
                for i in seq1:
                    x += i
                asa = x / ckLengt
                asal.append(asa)

            return asal

        # 判断一个fasta文件是否只有一条序列，返回一个布尔值
        def is_single_sequence(file_path):
            with open(file_path, "r") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            return len(records) == 1

        def calculate():
            self.ui.runButton.setEnabled(False)
            so.progress_update.emit(0)
            errNum = 3
            try:
                # 判断参考序列是否正确，只有一条或者文件格式正确
                filename1 = ''
                refilename = self.filePaths_ref[-1]
                if is_single_sequence(refilename):
                    filename1 = self.filePaths_alfa[-1]
                else:
                    errNum = 0
                    a1 > b1  # 故意不让后续代码运行

                so.progress_update.emit(1)

                if is_single_sequence(filename1):
                    errNum = 5
                    a1 > b1  # 故意不让后续代码运行
                else:
                    for record in SeqIO.parse(filename1, "fasta"):
                        if "X" in record.seq or "Z" in record.seq or "U" in record.seq or "B" in record.seq:
                            errNum = 4
                            a1 > b1  # 故意不让后续代码运行

                so.progress_update.emit(2)

                try:
                    pRl = [gswnPv(seq_record.seq) for seq_record in SeqIO.parse(filename1, "fasta")]
                except Exception:
                    errNum = 1
                so.progress_update.emit(3)

                # 使用列表推导式将每个子列表的行变为列，列变为行
                transposed = [[row[i] for row in pRl] for i in range(len(pRl[0]))]

                # pv值列表
                pvlst = [len(set(i)) for i in transposed]

                # print("PV：", pvlst)
                pv1 = re1pv(pvlst)  # 得到1/PV
                pvn = nor(pv1)  # 得到pv的归一化
                # print("输出1/PV归一化值：", pvn)

                # 存储pv和ASA
                data_pvASA = []

                line1 = [((str(i)) + ' to ' + (str(i + ckLengt - 1))) for i in range(1, len(pRl[0]) + 1)]
                line1.insert(0, ' ')

                data_pvASA.append(line1)
                # 想要创建一个新的、独立的列表，请使用切片操作符[:]来复制原始列表
                linePv = pvlst[:]
                linePv.insert(0, 'PV')
                data_pvASA.append(linePv)

                so.progress_update.emit(4)

                # 输入RSA
                try:
                    with open(self.filePaths_rsa[-1], 'r') as f:
                        data = f.read()

                    numbers = []
                    for line in data.split('\n'):
                        if not line.startswith('>'):
                            numbers += re.findall(r'\d+', line)

                    num = [int(x) for x in numbers]

                    asa = gswnAsa(num)
                    asa1 = asa[:]  # 想要创建一个新的、独立的列表，请使用切片操作符[:]来复制原始列表
                    asan = nor(asa)
                    # print("归一化ASA值：", asan)
                except Exception:
                    errNum = 2

                asa.insert(0, 'ASA')
                data_pvASA.append(asa)

                try:
                    filename = "./result/pv_ASA.xlsx"

                    if os.path.isfile(filename):
                        # 文件已经存在，重命名它
                        i = 1
                        while os.path.isfile(f'{filename}.{i}.xlsx'):
                            i += 1
                        os.rename(filename, f'{filename}.{i}.xlsx')

                    # 创建一个工作簿对象
                    workbook = openpyxl.Workbook()

                    # 获取活动工作表对象
                    worksheet = workbook.active

                    # 写入数据到单元格中
                    for row in data_pvASA:
                        worksheet.append(row)

                    # 保存工作簿
                    workbook.save('./result/pv_ASA.xlsx')
                    so.progress_update.emit(5)
                except Exception:
                    errNum = 6
                    a1 > b1

                # w = []
                wlst = list(np.add(pvn, asan))
                # print("输出W参数", wlst)

                # 列表按照值的大小进行降序排列，并记录降序之前的位置
                w_index = list(enumerate(wlst))
                w_sorted = sorted(w_index, key=lambda x: x[1], reverse=True)

                # print(w_sorted)

                # 计算序列长度，序列长度大于 1000 的前 100 名、序列长度在 500 和 1000 之间的前 75 名以及长度低于 500 的前 50 名
                seqlen = 0
                w_select = []  # （XX,XX）第一个是原索引，第二个是W数值
                for record in SeqIO.parse(filename1, "fasta"):
                    # print(len(record.seq))
                    seqlen = len(record.seq)
                    break
                if seqlen > 1000:
                    w_select = w_sorted[:100]
                elif 500 <= seqlen <= 1000:
                    w_select = w_sorted[:75]
                else:
                    w_select = w_sorted[:50]
                # print("降序排列，并选出：", w_select)

                # 弄一个存放降序排列的W值的文件
                wtrue = [[x[0] + 1, x[1]] for x in w_select]
                wtrue.insert(0, ['Starting position of the peptides', 'W value'])

                try:
                    filename = "./result/w_value.xlsx"

                    if os.path.isfile(filename):
                        # 文件已经存在，重命名它
                        i = 1
                        while os.path.isfile(f'{filename}.{i}.xlsx'):
                            i += 1
                        os.rename(filename, f'{filename}.{i}.xlsx')

                    # 创建一个工作簿对象
                    workbook = openpyxl.Workbook()

                    # 获取活动工作表对象
                    worksheet = workbook.active

                    # 写入数据到单元格中
                    for row in wtrue:
                        worksheet.append(row)

                    # 保存工作簿
                    workbook.save('./result/w_value.xlsx')

                    so.progress_update.emit(6)

                except Exception:
                    errNum = 7
                    a1 > b1

                winl = [w[0] for w in w_select]  # w存放了筛选出W参数的数列位置
                # print(winl)
                win = sorted(winl, reverse=False)
                aa_position = [i + 1 for i in win]  # 得到真实的位置

                # print("氨基酸的起始位置:", aa_position)

                # e = |f -f0|   不能为 0，因为两个 12 长的肽不能出现在同一位置

                # 这个列表的分组标准是：如果当前数与前一个数的差大于等于ckLengt，则将当前数单独分为一组；否则将当前数加入到上一个组中
                def group_list(lst):
                    groups = []
                    group = [lst[0]]
                    for i in range(1, len(lst)):
                        if lst[i] - lst[i - 1] >= ckLengt:
                            groups.append(group)
                            group = [lst[i]]
                        else:
                            group.append(lst[i])
                    groups.append(group)
                    return groups

                aa_group = group_list(aa_position)
                # print("分组完成：", aa_group)
                h = [lst[-1] - lst[0] + ckLengt for lst in aa_group]  # h表示分组肽区的长度
                # print("分组肽区的长度", h)

                # 每个分组区域的 PV值与ASA值得到
                result = []  # 存放aa_group首位到末尾的全部区间值
                for tpl in aa_group:
                    result.append(list(range(tpl[0], tpl[-1] + 1)))

                pv1gl = [list(pv1[i - 1] for i in gr) for gr in result]  # 各分组肽区的1/PV
                asagl = [list(asa1[i - 1] for i in gr) for gr in result]  # 各分组肽区的ASA

                # 每个分组区域的 ASA 值定义为构成该区域的所有 12 长肽的 ASA 值的平均值。同样，将每个分组肽段的 PV 值定义为其组成的 12 个长度肽段的 PV 值的平均值
                pv1g_mean = [np.mean(x) for x in pv1gl]  # 各分组肽区的平均1/PV
                asag_mean = [np.mean(x) for x in asagl]  # 各分组肽区的平均ASA

                pv1gn = nor(pv1g_mean)
                # print("归一化1/PVg值：", pv1gn)
                asagn = nor(asag_mean)
                # print("归一化ASAg值：", asagn)
                hn = nor(h)
                # print("归一化hn值：", hn)
                so.progress_update.emit(7)

                alst = []
                blst = []
                clst = []
                slst = []
                for i in range(len(hn)):
                    a = math.sqrt(pv1gn[i] ** 2 + hn[i] ** 2 + pv1gn[i] * hn[i])
                    alst.append(a)

                    b = math.sqrt(hn[i] ** 2 + asagn[i] ** 2 + hn[i] * asagn[i])
                    blst.append(b)

                    c = math.sqrt(asagn[i] ** 2 + pv1gn[i] ** 2 + asagn[i] * pv1gn[i])
                    clst.append(c)

                    s = (a + b + c) / 2
                    slst.append(s)

                so.progress_update.emit(8)
                abc_area = []  # 三角形面积
                for k in range(len(slst)):
                    area = math.sqrt(slst[k] * (slst[k] - alst[k]) * (slst[k] - blst[k]) * (slst[k] - clst[k]))
                    abc_area.append(area)
                # print("ΔABC的面积", abc_area)

                area_sorted = sorted(enumerate(abc_area), key=lambda x: x[1], reverse=True)
                area_top_50percent = area_sorted[:int(len(area_sorted) * 0.5)]

                # print("ΔABC的面积前50%", area_top_50percent)

                abcindex = [lst[0] for lst in area_top_50percent]

                # print("排名前50%的", abcindex)

                abc_50index = [result[x] for x in abcindex]
                # 得到氨基酸真实的起始位置和结束位置
                aa_start_end = []

                for aa in abc_50index:
                    if len(aa) > 1:
                        aa[-1] = aa[-1] + ckLengt - 1
                        aa_start_end.append([aa[0], aa[-1]])
                    else:
                        aa_start_end.append([aa[0], aa[-1] + ckLengt - 1])
                # print("前50%的肽区开始和结束位置", aa_start_end)
                so.progress_update.emit(9)

                # 存储最终结果
                data_result = []

                ref = []
                refilename = self.filePaths_ref[-1]

                for seq_record in SeqIO.parse(refilename, "fasta"):
                    ref = seq_record.seq

                aresult = [ref[aa_position[0] - 1:aa_position[1]] for aa_position in aa_start_end]
                aa_result = [str(seq) for seq in aresult]

                aa_start_end_result = [','.join(map(str, i)) for i in aa_start_end]
                # aa_start_end_result.insert(0, '前50%的肽区开始和结束位置')
                data_result.append(aa_start_end_result)
                data_result.append(aa_result)

                try:
                    filename = "./result/result.xlsx"
                    if os.path.isfile(filename):
                        # 文件已经存在，重命名它
                        i = 1
                        while os.path.isfile(f'{filename}.{i}.xlsx'):
                            i += 1
                        os.rename(filename, f'{filename}.{i}.xlsx')

                    # 创建一个工作簿对象
                    workbook = openpyxl.Workbook()

                    # 获取活动工作表对象
                    worksheet = workbook.active

                    # 写入数据到单元格中
                    for row in data_result:
                        worksheet.append(row)

                    # 保存工作簿
                    workbook.save('./result/result.xlsx')

                    # self.ongoing = False
                    so.progress_update.emit(10)
                except Exception:
                    errNum = 8
                    a1 > b1
                self.ui.runButton.setEnabled(True)
                # 发出完成信号
                so.finished.emit(data_result)
                # aa_result.insert(0, 'Sequence')【；。 data_result.append(aa_result)

            except Exception:
                so.error_window.emit(errNum)
                self.ui.runButton.setEnabled(True)



        # 创建一个线程对象，指定要执行的函数
        worker = Thread(target=calculate)
        # 启动线程
        worker.start()
        # # 等待子线程结束
        # worker.join()


if __name__ == '__main__':
    # 创建一个实例化对象，QApplication 提供了整个图形界面程序的底层管理功能
    app = QApplication([])
    # 加载 icon
    app.setWindowIcon(QIcon('logo.png'))
    # 创建一个实例化对象，接受定义的所有参数配置
    window = MainWindow()
    # 展示窗口及其所有的控件
    window.ui.show()
    # 进入事件处理循环（没有此段代码，窗口将会一闪而过）
    app.exec_()
