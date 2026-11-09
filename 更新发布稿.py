"""Rebuild data, five figures and numerical Markdown sections from one configuration.

Run with Python 3.10+; requires matplotlib for figures.
The opening, review framing, extensions and numerical sections are maintained here.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import re
from pathlib import Path

from 复算模型 import Model, export_data, load_config, mass_levels, run_checks
from 绘制图表 import draw_figures


def table(header: list[str], rows: list[list]) -> str:
    lines = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    lines.extend('| ' + ' | '.join(map(str, row)) + ' |' for row in rows)
    return '\n'.join(lines)


REVIEW_OPENING = r'''# 骑行传动、功率与阻力的数学模型：理论综述与可复算分析

*Mathematical Models of Cycling Transmission, Power and Resistance: A Theoretical Review with Reproducible Analyses*

## 摘要

骑行速度由道路与环境造成的机械功率需求、传动系统的转速及扭矩关系，以及骑手可持续输出共同决定。本文以已有道路骑行功率模型为基础，统一讨论齿比、踏频、曲柄与轮端功率、坡度、空气阻力、滚动阻力和加速惯性之间的数学关系；比较参数辨识、传动效率和个体踏频研究的建模目标与适用边界。在一维道路假设下，对巡航、匀加速与恒功率加速进行可复算分析，并以不同总质量和阻力情景，以及单前盘、双前盘的档位集合给出数值示例。基础物理关系采用已有理论，复合式与数表为本文整理和复算；数值校验与实测验证分别讨论。

关键词：骑行动力学；传动比；机械功率；踏频；空气阻力；滚动阻力；参数辨识

## 研究范围与文献脉络

本文是一篇叙述性理论综述与可复算分析，目标是把分散在骑行动力学、传动效率和生理输出研究中的关系放在一致的符号与单位体系下，解释它们如何组合、何时需要更复杂的模型。文献按建模问题组织，并沿代表性论文的参考文献扩展；不将这次范围性检索称为穷尽文献的系统综述，也不把已发表的基本关系主张为首次提出。

Olds等的道路表现模型[18]与Martin等的实测验证[1]构成机械功率需求的主要基础。Bos等[19]已将坡度、功率、踏频、齿比与垂直爬升速度联立；Danek等[20]把功率计和GPS测量用于阻力参数估计。Reed等[21]、Dunst等[22]及Dusling[24]分别从个体踏频、疲劳与加速过程研究骑手输出和齿比的耦合。传动损失模型[12][23]则解释了为什么相同几何齿比的不同零件组合也可能有不同实际表现。

| 模型分支 | 代表文献 | 主要问题 | 本文的整理与应用 |
|---|---|---|---|
| 道路功率需求与验证 | Olds等1995[18]；Martin等1998[1] | 阻力、势能和动能如何对应机械功率 | 统一符号、推导简化式并复算；不继承其他实验的预测精度 |
| 气动参数与测试方法 | Grappe等1997[13]；Malizia与Blocken2020[25] | 姿态、装备及测试条件怎样改变气动结果 | 解释CdA与Cd的区别、参考值来源及测试可比性 |
| 爬坡、踏频与齿比 | Bos等2020/2021[19] | 给定功率时坡度、踏频、齿比与爬升速度的关系 | 作为传动几何与道路需求联立的已有先例 |
| 参数辨识与不确定性 | Danek等2021[20] | 如何从功率和速度测量估计阻力参数 | 讨论可辨识性、有效参数与独立验证；本文情景网格没有概率含义 |
| 个体生理与动态输出 | Reed等2016[21]；Dunst等2023[22]；Dusling2026[24] | 踏频、疲劳和加速如何影响可输出功率 | 对比外部需求与人体供给；本文数值程序尚未拟合个体生理模型 |
| 传动效率与齿片尺寸 | Dowd等2025[12]；Barnaby等2020[23] | 扭矩、转速与齿片尺寸怎样影响损失 | 解释固定η的适用边界；[23]为会议短摘要，[24]目前仅核对出版商摘要 |

本文将三个相互关联但不同的问题分开：**道路需要多少功率、传动将它映射为什么踏频和扭矩、骑手能否持续提供该输出。**单盘、双盘、多前盘或内变速可在相同运动学定义下描述；具体装置的可用档位、效率和换挡过程需另建约束。本文提供的一维数值程序计算前两类关系，生理供给、侧风、下坡制动和换挡瞬态作为扩展讨论。

**阅读顺序：**第二、三、五节给出传动、阻力和加速推导；第一、八节说明参数依据与不确定性；第四、六、七节提供质量分级、不同齿数组合与传动范围算例；计算附录给出积分过程和复算方法。

'''


GENERAL_TRANSMISSION = r'''### 4. 通用传动比与多前盘档位集合

更一般地，将传动比定义为驱动啮合时轮端角速度除以曲柄角速度：

\[
G=\frac{\omega_w}{\omega_c},\qquad
v=\frac{C}{2\pi}\omega_w,\qquad
V=0.06CnG.
\tag{G1}
\]

这些关系不要求只有一个前盘。普通外变速的第i个前盘和第j个后齿片给出G_ij=N_f,i/N_r,j；若在后齿片与轮端之间还有内变速倍率h_k=ω_w/ω_sprocket，则G_ijk=(N_f,i/N_r,j)h_k。无内变速时h=1。带传动也可按相应轮的节距齿数计算；其他机构可直接采用实测或厂家定义并换算后的有效G。齿比定义的方向必须一致，不能把相反定义的数据直接混用。

用A表示经过链线、机械兼容性及使用策略限制后允许的前后档位组合，则实际档位集合为：

\[
\mathcal G=\left\{\frac{N_{f,i}}{N_{r,j}}:(i,j)\in\mathcal A\right\}.
\tag{G2}
\]

单前盘是i只有一个取值的特例；双前盘增加了可选组合，但可能出现近似重复齿比。档位数量不等于互不重叠的速度区间数量。固定齿比、连续可变传动也可分别用单个G或可变G(t)描述其运动学，但其自由轮、效率和控制特性并不因此相同。

式(3)的扭矩关系采用忽略传动内部能量暂存的准稳态效率定义。实际链条传动可以进一步设η=η(i,j,n,τ_c,润滑状态)，同样的G并不意味着完全相同的损失[12][23]。换挡过程还包含啮合转换、输出波动及可能的内部惯性，不能只凭瞬时齿数比完整描述。

在固定η、CdA和骑手给定输出的理想比较中，50×25与40×20的G都等于2，因此同踏频对应相同速度，同功率和踏频对应相同曲柄扭矩。若现实中效率或可输出功率不同，系统表现也可能不同。

'''


GENERAL_DYNAMICS = r'''### 第五步：统一能量平衡与模型扩展

式(8)是稳态特例。对沿道路运动的一维人车系统，保持v>0并将阻力功率按规定符号记账，可写为：

\[
\eta P_c
=\frac{dK}{dt}+Mg\,v\sin\theta
+P_{rr}+P_a+P_{\rm other}+P_{\rm brake}.
\tag{E1}
\]

K为平动和所计入转动部分的动能；P_rr为滚阻损失；P_other可包含另行建模的轴承、车轮旋转气动阻力等；P_brake≥0为制动耗散。这里的η将正向曲柄输入折算到轮端，若需要单独研究回驱、自由轮或储能传动，必须重新明确相应能量流向。[1][18]

本文后续计算取P_other=P_brake=0，并采用固定η，因而得到第五节的巡航—加速式。更一般的式(E1)说明各分支模型如何衔接，而不表示当前程序已经验证了所有工况。将某个损失并入有效CdA或Crr后，不能再把同一损失单独加一次。

在一维有符号来流扩展中，令w沿行进反方向为正，则u=v+w，气动力的阻力方向分量可写为F_a=0.5ρCdA·u|u|，P_a=F_a v。当来流始终相对迎面时退化为式(7)；强顺风使u<0时，气流可能对系统做正功，不能继续使用无符号平方阻力。下坡时势能项也可能为负，计算出的负需求应由滑行或制动解释，不等同于骑手提供“负功率”。

侧风需要向量来流、偏航角及随偏航变化的气动系数；车轮旋转气动损失也可能要单独辨识[25]。这些扩展不在当前仅接受非负坡度与正面逆风的一维程序中计算。

### 第六步：把人体供给与外部需求相连接

外部需求模型给出Pc(v,s,w,...)，人体模型则描述P_available(n,T,z)，其中T为持续时间，z表示疲劳或其他个体状态。研究可持续表现时，两者需要满足Pc≤P_available，并通过式(G1)将n与v、G关联。个体踏频研究[21]、疲劳下动态踏频模型[22]及生理驱动齿比模型[24]为这类耦合提供不同路线；它们的实验对象与优化目标并不相同，不能直接拼成通用的“最佳踏频”。

本文保留给定功率情景进行力学比较，没有为读者拟合P_available，也不把平均曲柄扭矩直接解释成肌肉力、膝关节负荷或代谢成本。这样的分层可以让同一套力学关系用于单盘、双盘或其他传动实例，而不把装置差别与人体差别混为同一个参数。

'''


REVIEW_CONCLUSION = '''## 结论

已有道路骑行模型、传动运动学和能量守恒可以在统一符号下联立，用于解释巡航、加速以及不同齿比对应的踏频和曲柄扭矩。单盘与双盘使用相同的基础关系，差别体现在实际档位集合、范围、重叠、链线与效率等条件；具体齿数组合属于通用模型的应用实例。

力学功率需求、传动映射和人体可输出功率是三个需要分别建模再耦合的环节。固定输出功率下改变齿比不会凭空增加机械功率，但改变踏频可能影响真实骑手的可持续输出；理想效率相同也不代表实际零件损失完全相同。

本文的数表与图形为给定假设下的可复算结果。其用途是解释变量关系、比较情景和检查推导，引用和数值校验都不能替代参数辨识与独立道路验证。综合分析应始终保留单位、功率测量端点、参数来源、适用条件与不确定性含义。

'''


def parameters_section(m: Model, config: dict) -> str:
    grid = config['grid']
    text = r'''## 一、计算条件、参数依据与适用边界

本文讨论单人、无跟骑、无刹车、沿直线道路骑行；坡度与正面逆风在所考察路段内不变。只计算平路和上坡，不外推到下坡制动、强顺风或侧风。风速是沿道路反向的分量，不是气象站风力等级。

### 1. 质量只输入一个人车装备总质量

总质量M包含骑手、自行车、衣物、水、工具和行李，车轮实际质量也已经包含在内。稳态重力和滚阻计算只需要这个总量，不再拆分车重与人重。本版取消“W/kg骑手”列，所有W/kg均以人车装备总质量为分母。

@PARAMETERS@

**当前正文参考系统为@MASS@kg；质量分级从@MIN@到@MAX@kg，每@STEP@kg一档。**这只是查表网格，不是人群体重分布。超出范围或处于两档之间，可直接在JSON中填写准确总质量，或扩展分级范围。

同一速度、坡度、逆风、CdA和Crr下，功率对总质量是线性的，因此73kg的功率可用70、75kg结果插值：P(73)=0.4P(70)+0.6P(75)。固定功率反求速度则不是线性关系，应按准确质量重算；不宜直接线性插值速度。临时选用较高质量档能提高同一阻力假设下的功率估计，但不保证覆盖其他参数误差。

不同质量等级采用独立的CdA情景，**不假定重5kg就必然增加多少迎风面积**。总质量足以作为本模型的质量输入，却不足以确定体型、骑姿和风阻。

### 2. 风阻输入是CdA，不能直接用车架的Cd

Cd是无量纲阻力系数，A是参考迎风面积；模型直接使用两者乘积CdA，单位m²，包含骑手与整车。仅有裸车架Cd、未说明参考面积的数字，或某款车在某条件下“节省多少瓦”，均不能直接当作人车CdA。

旧稿0.32m²没有可追溯的定量选值依据，本版不再将它用作初始基准。Grappe等[13]在室外自行车场，对同一名骑手及传统公路车测得UP（手握把横上部）、DP（下把）、AP（肘部支撑在休息把垫上）姿态的CdA分别为0.299、0.276、0.262m²。姿态比较实际只使用一名骑手，不能据此推断全体骑手平均值。本文的初始风阻情景按以下步骤构造：

| 情景 | CdA | 选值过程 | 证据性质 |
|---|---:|---|---|
| 低风阻参考 | 0.26m² | AP实测0.262取两位小数 | 对一个实测结果取整；不保证普通弯把姿态可达到 |
| 基准参考 | 0.30m² | UP实测0.299取两位小数 | 对一个实测结果取整；不是现代公路车人群均值 |
| 高风阻情景 | 0.40m² | 原文称0.299比既往UP报告低23%，0.299/(1−0.23)≈0.388，再向上选0.40 | 本文的情景构造，不是新增实测、标准上限或置信区间 |

原文用近似匀速的功率数据拟合总阻力与速度平方的关系，得到气动项。其原始关系可写成RT=Rr+0.5ρCdA·v²：拟合斜率k后CdA=2k/ρ。本文采用论文公布的反演结果，没有原始逐点数据来重新估计其统计误差。这里说明的是“实测→辨识→取整→情景比较”的过程，**不能从体重和车架名称唯一推导出某个人的CdA**。若修改配置中的CdA，则其来源需自行补充，不自动继承上述实测依据。

### 3. 真实车型的测试条件与可比性

| 真实车型 | 可追溯官方资料 | 对本模型的用途与限制 |
|---|---|---|
| Specialized Tarmac SL8、SL7、Venge | SL8白皮书[16]第4–7、21页：45km/h相关偏航分布；验证比较控制轮组、轮胎和车把堆高等 | 支持比较同一协议下的系统差异；本文未取得可直接移植的、统一零偏航绝对人车CdA数值，不给车型分配0.26/0.30/0.40标签 |
| Trek Madone Gen 8、Madone Gen 7、Émonda | 官方气动说明[17]：采用踩踏假人；坡度模拟统一轮胎、轮组及滚阻 | 说明车型气动测试需计入骑手及配置；不同厂商协议不能直接拼成绝对CdA排行榜 |

若以后取得同一骑手、姿态、空气密度、零偏航条件下两配置的**气动功率差**，无风时可用ΔCdA=2ΔP_air/(ρv³)换算；若差值明确是曲柄端且η相同，则先用ΔP_air=ηΔP_crank。只提供完赛时间差、偏航加权差或未说明功率测量端点时，不能直接套用这个式子。

### 4. 滚阻系数的实测参照与基准选值

BRR[14]公开测试的Continental Grand Prix 5000 S TR 28黑边样本，在72psi/5.0bar、实测宽29mm、无内胎加20ml补胎液条件下，单胎42.5kg载荷、页面约29km/h，给出9.7W和Crr=0.00291。测试协议[15]进一步说明实际滚筒速度28.8km/h、滚筒直径77cm、表面菱形纹理、温度21.5–22.5°C。按Crr=P/(Nv)复核：9.7/[42.5×9.80665×(28.8/3.6)]≈0.00291，与公开值一致（功率显示值经过取整）。

这提供了**已知产品、胎压、载荷与设备条件下的可追溯参照**，但它不是该轮胎在任意沥青路面的实测常数。不能用一个未经验证的固定倍率把滚筒值转换成所有真实路面的值，也不应把5bar解释为所有质量等级的推荐胎压。

Martin等[1]印刷页282转述Kyle(1988)十种高压开口胎在光滑沥青上的Crr范围0.0027–0.0040，原论文选均值0.0032。本文保留0.004作为初始基准，是明确选择该历史报道范围的上端，**不是声称现代轮胎平均值等于0.004**，也不保证它对所有道路都足够保守。此处是通过Martin转引，未把它包装成本文直接核验过的Kyle原始试验数据。

| 滚阻情景 | Crr | 用法 |
|---|---:|---|
| 指定轮胎滚筒参考 | 0.00291 | BRR上述样本的公开测量值；用于产品条件参照 |
| 文献范围上端基准 | 0.004 | 明确选用历史光滑沥青报道范围上端；用于基准比较 |
| 较高滚阻敏感性 | 0.006 | 本文定义0.004×1.5，检验结论对滚阻增加的敏感性；不对应已验证产品或路面等级 |

模型中Crr作用于**两轮合计法向载荷Mg·cosθ**。若前后轮的有效Crr相同，不再把Crr乘2；若不同，则Crr_eq=(Crr_f·N_f+Crr_r·N_r)/(N_f+N_r)，这里N_f、N_r指前后轮载荷，勿与后文齿数符号混用。本文仍用单一有效Crr，忽略载荷、胎压和速度变化引起的系数变化。

### 5. 怎样理解结果的可信范围？

数据按“总质量×CdA×Crr×道路条件”分别计算，读者应先选质量，再比较适用阻力情景。九组阻力组合是条件分析，**不是95%置信区间**；它们没有已知抽样分布，也不保证覆盖所有真实工况。增加引用能够使参数可追溯，不能替代对本人骑姿与道路的实测校准。

若要辨识自己的参数，可在近似无风、平路、无加速、同一姿态的多个速度点记录稳定功率与速度，在已知η假设下拟合ηP/v=MgCrr+0.5ρCdA·v²。截距对应MgCrr，斜率对应0.5ρCdA。需要足够速度跨度、重复测量并控制微小坡度与风；应另留数据检验预测误差，不能只报告拟合优度。未经独立测量时，η、风、坡度与功率计偏差都会混入拟合系数；结果更适合称为该配置的有效参数。

本文所有功率均为曲柄端机械功率；速度展示用km/h，动力学内部用m/s。模型不单独计算花鼓轴承阻力、车轮旋转气动阻力矩、腿部内部功、颠簸额外损耗及随档位变化的效率。等效转动质量仍只用于加速惯性，不能加进重力载荷。

'''
    rows = [
        ['人车装备总质量 M', f'{m.mass:g}kg', '单一输入；含车轮实际质量'],
        ['有效轮周 C', f'{m.circumference_m:g}m', '初始2.146m来自Garmin 700×30C参考[10]；最好实测'],
        ['人车阻力面积 CdA', f'{m.cda_m2:g}m²', '初始0.30的来源与取整过程见下文[13]'],
        ['有效滚阻系数 Crr', f'{m.crr:g}', '初始0.004为文献报道范围上端[1]；另列产品参照[14]'],
        ['空气密度 ρ', f'{m.rho:g}kg/m³', '固定计算条件；不代表实时当地空气'],
        ['传动效率 η', f'{m.eta:g}', '固定近似；非人体效率[12]'],
        ['重力加速度 g', f'{m.g:g}m/s²', '计算常数'],
        ['转动等效附加质量 m_rot', f'{m.rotational_equivalent_kg:g}kg', '仅用于加速；未测轮组惯量'],
    ]
    return (text.replace('@PARAMETERS@', table(['参数', '当前配置', '性质与来源'], rows))
            .replace('@MASS@', f'{m.mass:g}').replace('@MIN@', str(grid['mass_min_kg']))
            .replace('@MAX@', str(grid['mass_max_kg'])).replace('@STEP@', str(grid['mass_step_kg'])))


def cruise_section(m: Model, config: dict) -> str:
    result = '## 四、巡航数值与总质量分级速查\n\n'
    result += f'除另有说明，下列正文数表使用总质量{m.mass:g}kg、CdA={m.cda_m2:g}m²、Crr={m.crr:g}及配置中的其余参数。\n\n'
    result += '### 1. 无风不同坡度的巡航功率\n\n单位W；由式(8)计算并取整。\n\n'
    result += table(['坡度', '10km/h', '15km/h', '20km/h', '25km/h', '30km/h', '35km/h'],
                    [[f'{pct}%', *[round(m.cruise(v, pct / 100)) for v in (10, 15, 20, 25, 30, 35)]] for pct in (0, 2, 4, 6, 8, 10)])
    aero, rolling = m.forces(35)[2] * (35 / 3.6) / m.eta, m.forces(35)[1] * (35 / 3.6) / m.eta
    result += f'\n\n平路无风35km/h：气动约{aero:.1f}W、滚阻约{rolling:.1f}W，总计约**{m.cruise(35):.1f}W**，均已折算到曲柄端。高功率格只表示力学需求，不表示骑手可持续达到。\n\n'
    result += '![无风不同坡度的匀速巡航功率曲线](图1_匀速巡航_坡度.png)\n\n### 2. 平路不同正面逆风的巡航功率\n\n'
    result += table(['地速', '无风', '10km/h正逆风', '20km/h正逆风'],
                    [[f'{v}km/h', *[f'{m.cruise(v, 0, w):.0f}W' for w in (0, 10, 20)]] for v in (20, 25, 30, 35, 40, 45)])
    result += '\n\n![平路正面逆风的匀速巡航功率曲线](图2_匀速巡航_逆风.png)\n\n### 3. 每5kg一档的总质量速查\n\n'
    result += f'这里固定CdA={m.cda_m2:g}m²、Crr={m.crr:g}、无风，单独观察质量变化。前三列为巡航所需曲柄功率，末两列为170W爬坡稳态速度。更换风阻、滚阻或功率请筛选随附CSV。\n\n'
    result += table(['总质量', '平路25km/h', '平路30km/h', '平路35km/h', '6%坡170W', '8%坡170W'],
                    [[f'{mass}kg', *[f'{replace(m, total_mass_kg=mass).cruise(v):.0f}W' for v in (25, 30, 35)],
                      *[f'{replace(m, total_mass_kg=mass).speed_at_power(170, grade):.2f}km/h' for grade in (.06, .08)]] for mass in mass_levels(config)])
    result += '\n\n![总质量与风阻情景对比](图4_质量与风阻分级.png)\n\n### 4. W/kg只按人车装备总质量归一化\n\n'
    result += r'''将式(8)除以实际总质量M：

\[
\boxed{\frac{P_{\mathrm{cruise}}}{M}=\frac v\eta\left[g\frac{s+C_{rr}}{\sqrt{1+s^2}}+\frac12\rho\frac{C_dA}{M}(v+w)^2\right]}
\tag{10}
\]

来源与推导：式(8)按总质量归一化[1][6][7]。M不包含额外转动等效质量。公式仍保留CdA/M，因此不同质量、体型和骑姿不能共用一张无条件的W/kg表。

'''
    power = m.cruise(15, .06)
    return result + f'本配置6%坡15km/h约需{power:.1f}W，即{power/m.mass:.2f}W/kg总质量。本版未提供以骑手体重为分母的训练指标。\n\n'


def acceleration_numbers(m: Model) -> str:
    result = '### 4. 平路无风，从20加到35km/h\n\n'
    result += f'当前总质量{m.mass:g}kg，M_eff={m.effective_mass:g}kg。以下为匀加速，不能当作恒功率过程。\n\n'
    results = [(t, m.constant_accel(20, 35, t)) for t in (10, 15, 30)]
    result += table(['用时', '加速度', '过程平均曲柄功率', '到35前一刻仍加速', '到35后巡航'],
                    [[f'{t}秒', f'{a["acceleration_m_s2"]:.3f}m/s²', f'{a["mean_total_power_w"]:.1f}W',
                      f'{a["end_power_w"]:.1f}W', f'{a["final_cruise_power_w"]:.1f}W'] for t, a in results])
    a = results[0][1]
    result += f'\n\n三种过程的动能增量都约{a["kinetic_energy_increase_j"]/1000:.3f}kJ，额外曲柄做功约{a["kinetic_energy_increase_j"]/m.eta/1000:.3f}kJ；平均阻力功率都约{a["mean_resistance_power_w"]:.1f}W，因为它们在相同速度区间线性提速。\n\n'
    result += '平均功率低于终点巡航需求并不矛盾：过程中功率逐渐增加，终点前仍需额外加速功率。慢加速会走更远、克服阻力更久，不能把平均功率较低等同于全程机械做功较少。\n\n![加速结束后转入巡航](图3_加速与巡航.png)\n\n### 5. 恒定功率时，需要多久才能到35km/h？\n\n'
    rows = []
    for power in (150, 170, 201, 250, 300):
        reached = power > m.cruise(35)
        rows.append([f'{power}W', f'{m.speed_at_power(power):.2f}km/h',
                     f'{m.constant_power_time(20, 35, power):.2f}秒' if reached else '不能有限时间达到'])
    result += table(['恒定曲柄功率', '平路无风稳态速度', '20→35km/h用时'], rows)
    return result + '\n\n依据式(12)分离变量积分，详见计算附录。目标需求高于给定功率时不可达；恰好相等时渐近接近，不报告伪造的有限到达时间。未模拟换挡中断和人体输出波动。\n\n'


def gearing_section(m: Model, config: dict) -> str:
    rear = max(config['grid']['cassette_teeth'])
    result = '## 六、把功率映射为齿比、踏频和扭矩\n\n给定道路条件和曲柄功率，先由式(8)解稳态速度，再用式(1)、(2)算踏频和扭矩。170W只是示例，不代表FTP或任意踏频下的可持续输出。\n\n'
    rows = []
    for grade in (.06, .08):
        speed = m.speed_at_power(170, grade)
        for front in (50, 46, 44):
            rpm = m.cadence(front, rear, speed)
            rows.append([f'{grade:.0%}坡、无风、170W', f'{front}×{rear}', f'{speed:.2f}km/h', f'{rpm:.1f}rpm', f'{m.crank_torque(170, rpm):.1f}N·m'])
    result += table(['场景', '档位', '稳态速度', '踏频', '曲柄平均扭矩'], rows)
    result += '\n\n同坡、同功率、同阻力假设下，齿比改变踏频和扭矩，稳态速度相同。上表每一行表示一个独立的前后齿数组合，可出现在单盘或多盘系统中；小前盘改善轻档，并不保证爬坡时达到个人设定踏频。\n\n### 单前盘与双前盘的几何算例\n\n'
    rear_min = min(config['grid']['cassette_teeth'])
    speed = m.speed_at_power(170, .08)
    examples = []
    for name, fronts in [('44T单前盘', [44]), ('50/34T双前盘', [50, 34]), ('52/36T双前盘', [52, 36])]:
        smallest = min(fronts) / rear
        largest = max(fronts) / rear_min
        examples.append([name, f'{smallest:.3f}', f'{largest:.3f}', f'{100*largest/smallest:.1f}%',
                         f'{m.cadence(min(fronts), rear, speed):.1f}rpm',
                         f'{m.speed_at_cadence(max(fronts), rear_min, 90):.1f}km/h'])
    result += table(['示例', '最小G', '最大G', '理论跨度', '8%坡170W最低档踏频', '最高档90rpm速度'], examples)
    result += f'\n\n三个例子使用相同{rear_min}–{rear}后齿片范围、总质量和阻力参数。双前盘例子在几何上同时提供更低的最低档和更高的最高档；表内速度仍受实际可用功率约束。这不是产品装车建议，具体允许组合、重复齿比、链线和后拨容量还要另查。当前附带的批量CSV以配置中的前盘候选逐个列行，不将它们全部视为同一辆车上的前盘。\n\n### 功率下降后的轻档余量\n\n'
    result += table(['8%坡无风功率', '稳态速度', f'44×{rear}踏频'],
                    [[f'{p}W', f'{m.speed_at_power(p, .08):.2f}km/h', f'{m.cadence(44, rear, m.speed_at_power(p, .08)):.1f}rpm'] for p in config['grid']['powers_w']])
    result += '\n\n请把功率与持续时间联系起来，例如按自己的实测10分钟、30分钟或疲劳后可持续功率重新计算；这里150/170/210W是数值情景，并未被指定为某种生理能力等级。保底爬坡档应考虑长坡后段和疲劳余量。70rpm也只是可修改的个人约束。\n\n### 高速端的纯几何上限\n\n'
    result += table(['前盘×最小齿片', '90rpm对应速度'],
                    [[f'{front}×{min(config["grid"]["cassette_teeth"])}', f'{m.speed_at_cadence(front, min(config["grid"]["cassette_teeth"]), 90):.1f}km/h'] for front in (50, 48, 46, 44, 42)])
    return result + '\n\n这只是踏频与速度的对应关系，不是最高可达到或持续的速度，更不表示给定功率能在平路维持该速度。\n\n'


def range_section(m: Model, config: dict) -> str:
    grid = config['grid']
    low, high = grid['cadence_min_rpm'], grid['cadence_max_rpm']
    rear_min, rear_max = min(grid['cassette_teeth']), max(grid['cassette_teeth'])
    slow = m.speed_at_power(170, .08)
    front_min = 45 * rear_min / (.06 * m.circumference_m * high)
    front_max = slow * rear_max / (.06 * m.circumference_m * low)
    required_span = 45 * low / (slow * high)
    result = r'''## 七、离散传动的范围、档位间隔与可行性

### 1. 任意传动系统的踏频可行区间

希望给定速度V时踏频在n_min至n_max之间，由式(1)移项得到：

\[
\boxed{\frac{V}{0.06Cn_{\max}}\le G\le\frac{V}{0.06Cn_{\min}}}
\tag{15}
\]

实际允许的传动比集合𝒢中存在某个档位落入区间，才满足该工况的踏频要求。单前盘检查各个后齿片，双前盘检查允许的前后组合；满足转速仍不保证生理上舒适或可持续输出。

'''
    result += f'### 2. 固定一个前盘的边界算例\n\n用当前参考系统举例：45km/h时不超过{high:g}rpm，同时8%坡170W时不低于{low:g}rpm；此算例仅允许一个前盘，后齿片限定{rear_min}–{rear_max}。爬坡速度为{slow:.4f}km/h，反推前盘：\n\n'
    result += '\\[\nN_f\\ge\\frac{' + f'45\\times{rear_min}' + '}{' + f'0.06\\times{m.circumference_m:g}\\times{high:g}' + '}\\approx' + f'{front_min:.2f}' + '\n\\]\n\n'
    result += '\\[\nN_f\\le\\frac{' + f'{slow:.4f}\\times{rear_max}' + '}{' + f'0.06\\times{m.circumference_m:g}\\times{low:g}' + '}\\approx' + f'{front_max:.2f}' + '\n\\tag{16}\n\\]\n\n'
    result += ('两端没有交集，当前飞轮无法同时满足这两个条件。' if front_min > front_max else '两端存在交集，但还要检查可选整数前盘和中间实际档位。') + '\n\n'
    result += r'''同一前盘下，飞轮决定跨度：

\[
\frac{G_{\max}}{G_{\min}}=\frac{N_f/N_{r,\min}}{N_f/N_{r,\max}}=\frac{N_{r,\max}}{N_{r,\min}}
\tag{17}
\]

进一步将高速端最低齿比与爬坡端最高齿比相除，可得连续齿比范围的必要条件：

\[
\frac{G_{\max}}{G_{\min}}\ge\frac{V_{\rm high}n_{\min}}{V_{\rm climb}n_{\max}}
\]

'''
    result += f'本例至少需要{required_span*100:.1f}%跨度，单前盘配{rear_min}–{rear_max}只有{rear_max/rear_min*100:.1f}%。在保持单前盘结构时，替换前盘只会整体移动范围，不会扩大跨度；增加第二个前盘则改变可用齿比集合。这个判据并不检查中间档位间隔和可安装零件。\n\n'
    result += r'''### 3. 多前盘的范围与重叠

如果前后所有组合均允许，且无内变速，则：

\[
\frac{G_{\max}}{G_{\min}}
=\frac{N_{f,\max}}{N_{f,\min}}
\frac{N_{r,\max}}{N_{r,\min}}.
\tag{G3}
\]

若某些组合因链线或兼容性被排除，应直接在式(G2)允许集合中取最大值与最小值，不能无条件使用式(G3)。例如50/34前盘配11–34后齿片的几何跨度为(50/34)×(34/11)≈454.5%，这一结果以两个极端组合均允许为前提。

将允许齿比去重后排序，可以观察不同前盘速度区间的重叠及相邻比值间隔；实际换挡路径还取决于前后拨操作顺序，按齿比排序相邻的两个组合并不一定能通过一次后拨操作互换。

### 4. 常用速度下的具体档位

'''
    result += f'配置中的飞轮齿数示例为：**{", ".join(map(str, grid["cassette_teeth"]))}T**。不是只有首尾两片；完整组合可在JSON修改。35km/h下各档踏频如下：\n\n![常用巡航速度的档位踏频](图5_常用巡航档位.png)\n\n'
    result += table(['前盘', '后齿片', '35km/h踏频', f'位于{low:g}–{high:g}rpm'],
                    [[str(front), str(rear), f'{m.cadence(front, rear, 35):.1f}', '是' if low <= m.cadence(front, rear, 35) <= high else '否']
                     for front in (42, 44, 46, 50) for rear in grid['cassette_teeth'] if 12 <= rear <= 19])
    result += '\n\n车速不变时，从后齿片r₁换到更大r₂，踏频相对原值增加r₂/r₁−1。例如90rpm时15→17T上升13.3%，变成102rpm。若在固定踏频下比较速度，百分比则是另一个分母，不能混用。CSV同时给出各档在设定踏频上下限内的速度区间。\n\n'
    return result + '数学可行不等于装车兼容。R8150官方表列出11–30和11–34，不含11–36[11]；盘片齿形、链线、保链、后拨总容量与塔基都应另核实。双前盘还需检查前拨匹配及齿差容量，单前盘也有相应保链要求。“很少使用”的保底爬坡档不自动等于浪费。\n\n'


def limitations_section(m: Model, config: dict) -> str:
    rows = []
    for aero in config['aero_profiles']:
        rows.append([f'{aero["value"]:g}', *[f'{replace(m, cda_m2=aero["value"], crr=rolling["value"]).cruise(35):.1f}' for rolling in config['rolling_profiles']]])
    return ('## 八、参数敏感性与结果边界\n\n'
            f'以下为总质量{m.mass:g}kg、平路无风35km/h的曲柄功率（W），列分别对应三种Crr，所有其他参数相同。\n\n'
            + table(['CdA（m²）', *[f'Crr={p["value"]:g}' for p in config['rolling_profiles']]], rows)
            + '\n\n同一质量下也可能因姿态、轮胎和路况而产生明显差异。表内最小值与最大值只是所选组合的情景范围，没有概率含义，不可标作统计置信区间。不要用一个假定CdA对不同车架做实测性能排名。\n\n'
            + '模型能估算道路功率需求、条件稳态速度、踏频和扭矩以及传动范围边界；不能判断某个人能坚持多久、哪个踏频代谢最省、膝关节负荷是否适宜或跟骑收益。实际效率会随工作状态变化[12]，本版仍使用固定η。\n\n'
            + '稳态公式不能直接代入整趟平均速度、平均坡度。路线分析应按时间或路段分别计算，再积分能量；刹车和反复加速也要单独处理。\n\n')


def appendix(m: Model) -> str:
    return r'''# 计算附录：加速过程的平均功率与恒功率积分

本附录与主文使用同一配置。下列解析式基于固定道路条件、固定CdA、Crr、η以及匀加速假设，不是人体输出模型。

## 匀加速平均功率

正文式(8)写成：

\[
P_{\rm cruise}(v)=\frac{Av+Bv(v+w)^2}{\eta},\quad
A=Mg\frac{s+C_{rr}}{\sqrt{1+s^2}},\quad B=\frac12\rho C_dA.
\tag{A1}
\]

由v(t)=v₀+at，aT=v₁−v₀，对于v₁>v₀：

\[
\overline{v^k}=\frac1{v_1-v_0}\int_{v_0}^{v_1}v^k\,dv
=\frac{v_1^{k+1}-v_0^{k+1}}{(k+1)(v_1-v_0)}.
\tag{A2}
\]

\[
\overline v=\frac{v_0+v_1}{2},\quad
\overline{v^2}=\frac{v_0^2+v_0v_1+v_1^2}{3},\quad
\overline{v^3}=\frac{v_0^3+v_0^2v_1+v_0v_1^2+v_1^3}{4}.
\tag{A3}
\]

展开式(A1)并按时间平均：

\[
\overline P_{\rm resistance}=\frac{A\overline v+B(\overline{v^3}+2w\overline{v^2}+w^2\overline v)}{\eta}.
\tag{A4}
\]

加上动能增量得到：

\[
\boxed{\overline P_c=\overline P_{\rm resistance}
+\frac{(M+m_{\rm rot})(v_1^2-v_0^2)}{2\eta T}}
\tag{A5}
\]

基础关系与来源：主文[1][4][7][8][9]。相同速度区间、相同条件的匀加速过程有相同平均阻力功率，但不同T对应不同位移、阻力做功与额外平均加速功率。v₁=v₀时直接用巡航式，不能在(A2)里除以零。

## 恒功率加速

由正文式(12)移项（v>0）：

\[
\frac{dv}{dt}=\frac{\eta P_c/v-F_{\rm resistance}(v)}{M+m_{\rm rot}}.
\tag{A6}
\]

分离变量可计算从v₀到v₁的用时：

\[
T=\int_{v_0}^{v_1}\frac{(M+m_{\rm rot})v}{\eta[P_c-P_{\rm cruise}(v)]}\,dv.
\tag{A7}
\]

这里v₀>0、v₁>v₀，并要求整个速度区间内P_c>P_cruise(v)。在本文非负坡度/逆风及正阻力条件下，P_cruise随速度单调增加，可检查终点。若P_c低于终点需求则无法到达；恰好等于时渐近接近。数值代码以逐步加密的复合Simpson积分求用时，并检查收敛；太接近渐近边界而未收敛时明确报错，不输出虚假的到达时间。

从静止出发不能直接代入(A6)，需另引入有限曲柄扭矩、牵引力和控制过程。换挡导致的踏频、功率波动也未在此模拟。

## 当前配置的数值核对

''' + f'总质量{m.mass:g}kg、CdA={m.cda_m2:g}m²、Crr={m.crr:g}、η={m.eta:g}，平路无风20→35km/h：\n\n' + table(
        ['用时', '动能增量J', '平均阻力W', '平均总功率W', '末端功率W'],
        [[str(t), *[f'{m.constant_accel(20,35,t)[key]:.6f}' for key in ('kinetic_energy_increase_j', 'mean_resistance_power_w', 'mean_total_power_w', 'end_power_w')]] for t in (10,15,30)]) + f'\n\n35km/h巡航为{m.cruise(35):.6f}W。小数位仅用于复算一致性，不表示实测精度。\n\n## 校验与再生成\n\n`复算模型.py --check`只校验配置、功率反求速度及解析平均功率与梯形积分的一致性，不写文件。`test_复算模型.py`另以独立受力表达式、Simpson积分及RK4动态积分检查回归、单位和边界，包含无效输入、不可达恒功率目标与配置生效测试。\n\n`更新发布稿.py`从JSON生成全部CSV、五张图和正文数值章节。先校验，再生成；这些检查不等于实验验证，参数辨识与独立道路验证仍需实测。\n'


REFERENCES = '''

[13] Grappe F, Candau R, Belli A, Rouillon JD. *Aerodynamic drag in field cycling with special reference to the Obree’s position*. Ergonomics. 1997;40(12):1299–1311. DOI: 10.1080/001401397187388。摘要、方法§2.1–2.2与结果：姿态辨识使用单名骑手；UP/DP/AP的CdA为0.299/0.276/0.262m²。https://doi.org/10.1080/001401397187388 ；作者公开全文：https://www.researchgate.net/profile/Robin-Candau/publication/283892680_Grappe_1997_aerodynamic_resistance_in_field_cycling_with_special_reference_to_the_Obrees_position/links/5649969108ae54697fbe8361/Grappe-1997-aerodynamic-resistance-in-field-cycling-with-special-reference-to-the-Obrees-position.pdf

[14] Bicycle Rolling Resistance, Jarno Bierman. *Continental Grand Prix 5000 S TR 28*. 2023-01-06。公开Rolling Resistance Test Results表：72psi、单胎42.5kg载荷、约29km/h、9.7W、Crr=0.00291。https://www.bicyclerollingresistance.com/road-bike-reviews/continental-grand-prix-5000-s-tr-28

[15] Bicycle Rolling Resistance. *Test Info / FAQ*, §1.1–1.2。滚筒、载荷、温度与单胎/双胎功率说明。https://www.bicyclerollingresistance.com/the-test

[16] Specialized. *Tarmac SL8 White Paper*. 2023，第4–7、21页。测试速度与偏航分布、SL8/SL7/Venge比较控制条件。https://specialized.a.bigcontent.io/v1/static/ROAD-5793-White-Paper-2023-Tarmac-SL8

[17] Trek. *Madone Gen 8 aerodynamics*. “Weight vs aerodynamics”部分：踩踏假人及统一轮胎、轮组与滚阻的对比条件。https://www.trekbikes.com/in/en_IN/madone-aerodynamics/

[18] Olds TS, Norton KI, Lowe EL, Olive S, Reay F, Ly S. *Modeling road-cycling performance*. Journal of Applied Physiology. 1995;78(4):1596–1611. DOI: 10.1152/jappl.1995.78.4.1596。道路骑行功率需求的第一性原理模型。https://pubmed.ncbi.nlm.nih.gov/7615475/

[19] Bos L, Slawinski MA, Stanoev T. *On maximizing VAM for a given power: Slope, cadence, force and gear-ratio considerations*. arXiv:2006.15816，2020首投，2021年v5修订。参见式(1)、(14)及§3.3。https://arxiv.org/html/2006.15816v5

[20] Danek T, Slawinski MA, Stanoev T. *On modelling bicycle power-meter measurements*. arXiv:2103.09806，2021。讨论功率计/GPS数据、阻力参数反演及不确定性。https://arxiv.org/abs/2103.09806

[21] Reed R, Scarf P, Jobson SA, Passfield L. *Determining optimal cadence for an individual road cyclist from field data*. European Journal of Sport Science. 2016;16(8):903–911. DOI: 10.1080/17461391.2016.1146336。https://pmc.ncbi.nlm.nih.gov/articles/PMC4989856/

[22] Dunst AK, Hesse C, Ueberschär O. *The Concept of Optimal Dynamic Pedalling Rate and Its Application to Power Output and Fatigue in Track Cycling Sprinters—A Case Study*. Sports. 2023;11(1):19. DOI: 10.3390/sports11010019。https://pmc.ncbi.nlm.nih.gov/articles/PMC9864339/

[23] Barnaby GC, Yon J, Burgess S. *Sprocket Size Optimisation for Derailleur Racing Bicycles*. Journal of Science and Cycling. 2020;9(2):36。会议特刊短摘要，讨论齿片尺寸、损失与质量/惯量的权衡；不能据此复原完整方法。https://www.jsc-journal.com/index.php/JSC/article/view/581

[24] Dusling K. *A physiologically motivated analytic model for optimal bicycle gearing*. Proceedings of the Institution of Mechanical Engineers, Part P: Journal of Sports Engineering and Technology. OnlineFirst，2026-10-26. DOI: 10.1177/17543371261499757。本文目前仅依据出版商摘要概述其生理功率模型、运动解析式与加速齿比优化范围，未复现全文方法。https://journals.sagepub.com/doi/abs/10.1177/17543371261499757

[25] Malizia F, Blocken B. *Bicycle aerodynamics: History, state-of-the-art and future perspectives*. Journal of Wind Engineering and Industrial Aerodynamics. 2020;200:104134. DOI: 10.1016/j.jweia.2020.104134。综述自行车气动研究及测试背景，具体参数仍需与相应原始试验条件对应。https://doi.org/10.1016/j.jweia.2020.104134
'''


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    first = text.index(start)
    last = text.index(end, first + len(start))
    return text[:first] + replacement + text[last:]


def update_markdown(m: Model, config: dict, root: Path) -> None:
    path = root / '初稿.md'
    text = path.read_text(encoding='utf-8')
    text = REVIEW_OPENING + text[text.index('## 一、'):]
    transmission = text[text.index('## 二、'):text.index('## 三、')]
    transmission = transmission.split('### 4. 通用传动比与多前盘档位集合', 1)[0]
    transmission = transmission.replace('## 二、传动：齿比不是功率档位', '## 二、传动运动学与扭矩映射')
    text = replace_section(text, '## 二、', '## 三、', transmission.rstrip() + '\n\n' + GENERAL_TRANSMISSION)
    dynamics = text[text.index('## 三、'):text.index('## 四、')]
    dynamics = dynamics.split('### 第五步：统一能量平衡与模型扩展', 1)[0]
    text = replace_section(text, '## 三、', '## 四、', dynamics.rstrip() + '\n\n' + GENERAL_DYNAMICS)
    text = replace_section(text, '## 一、', '## 二、', parameters_section(m, config))
    text = replace_section(text, '## 四、', '## 五、', cruise_section(m, config))
    text = replace_section(text, '### 4. 平路无风，从20加到35km/h', '## 六、', acceleration_numbers(m))
    text = replace_section(text, '## 六、', '## 七、', gearing_section(m, config))
    text = replace_section(text, '## 七、', '## 八、', range_section(m, config))
    text = replace_section(text, '## 八、', '## 结论', limitations_section(m, config))
    text = replace_section(text, '## 结论', '## 参考文献与公式出处', REVIEW_CONCLUSION)
    # Keep sources [1]–[12] and their numbering, replacing only this script's additions.
    text = re.split(r'\n\[13\] ', text, maxsplit=1)[0].rstrip() + REFERENCES + '\n'
    path.write_text(text, encoding='utf-8')
    (root / '计算附录.md').write_text(appendix(m), encoding='utf-8')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('模型参数.json'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    m, config = load_config(args.config)
    run_checks(m)
    # Check the optional renderer before writing outputs.
    import matplotlib  # noqa: F401
    export_data(m, config, root)
    draw_figures(m, config, root)
    update_markdown(m, config, root)
    print('Publication data, five figures and Markdown numerical sections updated.')


if __name__ == '__main__':
    main()
