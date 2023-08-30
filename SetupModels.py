#execfile("FoldersName.py")
import brian2 as b2
import numpy as np
from FoldersName import *
from NeuronModelsVarDrive import *

Models = {0:'IzhiTypeI', 1:'IzhiTypeII', 2:'LIFModel', 3:'AdExModel', 4:'GuillemModel', 5:'IzhiTypeIIModified',6:'WangBuzsaki'}
SubModelsAdEx = {0:'TypeI', 1:'TypeIISaddleNode', 2:'TypeIIAndronovHopf'}

FoldersI0 = ['I0_LIF','I0_IzhiTyII','I0_IzhiTyI','I0_Guillem','I0_AdExTyI','I0_AdExTyIISN','I0_AdExTyIIAH','I0_IzhiTyIIMod']

FoldersCSV = ['CSV_LIF','CSV_IzhiTyII','CSV_IzhiTyI','CSV_Guillem','CSV_AdExTyI','CSV_AdExTyIISN','CSV_AdExTyIIAH','CSV_IzhiTyIIMod']

FoldersCSVTheta = ['CSVTheta_LIF','CSVTheta_IzhiTyII','CSVTheta_IzhiTyI','CSVTheta_Guillem','CSVTheta_AdExTyI','CSVTheta_AdExTyIISN','CSVTheta_AdExTyIIAH','CSVTheta_IzhiTyIIMod']

FoldersFigures = ['Fig_LIF','Fig_IzhiTyII','Fig_IzhiTyI','Fig_Guillem','Fig_AdExTyI','Fig_AdExTyIISN','Fig_AdExTyIIAH','Fig_IzhiTyIIMod']

FoldersFiguresTheta = ['FigTheta_LIF','FigTheta_IzhiTyII','FigTheta_IzhiTyI','FigTheta_Guillem','FigTheta_AdExTyI','FigTheta_AdExTyIISN','FigTheta_AdExTyIIAH','FigTheta_IzhiTyIIMod']

IextUnits = 1.*b2.mV

# Izhikevich Type I calibrated to CA1 FS PV+ neurons (from Ferguson et al. 2013)
if ChosenModel==0:
    Model = IzhiTypeI
    threshold_cond = threshold_condIzhiTypeI
    reset_cond = reset_condIzhiTypeI
    v_th = v_thIzTyI
    v_r = v_rIzTyI
    Cm = CmIzTyI
    uUnits = uUnitsIzTyI
    v0 = v0IzTyI
    u0 = u0IzTyI

    I0Folder = I0FolderIzhiTypeI
    if ThetaInput:
        CSVFolder = CSVFolderIzhiTypeITheta
        FiguresFolder = FiguresFolderIzhiTypeITheta
    elif VariableDrive:
        CSVFolder = CSVFolderIzhiTypeIVarDrive
        FiguresFolder = FiguresFolderIzhiTypeIVarDrive		
    else:
        CSVFolder = CSVFolderIzhiTypeI 
        FiguresFolder = FiguresFolderIzhiTypeI

# Izhikevich Type II resonator
elif ChosenModel==1:
    Model = IzhiTypeII 
    threshold_cond = threshold_condIzhiTypeII
    reset_cond = reset_condIzhiTypeII
    v_th = v_thIzTyII
    v_r = v_rIzTyII
    aIzTyII = aIzTyII
    uUnits = uUnitsIzTyII
    v0 = v0IzTyII
    u0 = u0IzTyII

    I0Folder = I0FolderIzhiTypeII
    if ThetaInput:
        CSVFolder = CSVFolderIzhiTypeIITheta
        FiguresFolder = FiguresFolderIzhiTypeIITheta
    elif VariableDrive:
        CSVFolder = CSVFolderIzhiTypeIIVarDrive
        FiguresFolder = FiguresFolderIzhiTypeIIVarDrive	
    else:
        CSVFolder = CSVFolderIzhiTypeII
        FiguresFolder = FiguresFolderIzhiTypeII

# Leaky Integrate and Fire
elif ChosenModel==2:
    Model = LIFModel
    threshold_cond = threshold_condLIF
    reset_cond = reset_condLIF
    v_th = v_thLIF
    v_r = v_rLIF
    v0 = v0LIF
    u0 = 1.# Dummy variable for LIF model

    I0Folder = I0FolderLIF
    if ThetaInput:
        CSVFolder = CSVFolderLIFTheta
        FiguresFolder = FiguresFolderLIFTheta
    elif VariableDrive:
        CSVFolder = CSVFolderLIFVarDrive
        FiguresFolder = FiguresFolderLIFVarDrive	
    else:
        CSVFolder = CSVFolderLIF
        FiguresFolder = FiguresFolderLIF

# Adaptive Exponential 
elif ChosenModel==3:
    Model = AdExModel
    threshold_cond = threshold_condAdEx
    reset_cond = reset_condAdEx
    wUnits = wUnitsAdEx
    v0 = v0AdEx
    u0 = w0AdExSN

    # AdEx Type I with Saddle Node bifurcation
    if ChosenSubmodelAdEx==0:
        glAdEx, CAdEx, ElAdEx, DeltaTAdEx, V_TAdEx, tau_wAdEx, aAdEx, bAdEx = glAdExSN, CAdExSN, ElAdExSN, DeltaTAdExSN, V_TAdExSN, tau_wAdExSN, aAdExSN, bAdExSN
        v_rAdEx = v_rAdExSNTypeI
        I0Folder = I0FolderAdExTypeI
        if ThetaInput:
            CSVFolder = CSVFolderAdExTypeITheta
            FiguresFolder = FiguresFolderAdExTypeITheta
        elif VariableDrive:
            CSVFolder = CSVFolderAdExTypeIVarDrive
            FiguresFolder = FiguresFolderAdExTypeIVarDrive	
        else:
            CSVFolder = CSVFolderAdExTypeI 
            FiguresFolder = FiguresFolderAdExTypeI
    # AdEx Type II with Saddle Node bifurcation
    elif ChosenSubmodelAdEx==1:
        glAdEx, CAdEx, ElAdEx, DeltaTAdEx, V_TAdEx, tau_wAdEx, aAdEx, bAdEx = glAdExSN, CAdExSN, ElAdExSN, DeltaTAdExSN, V_TAdExSN, tau_wAdExSN, aAdExSN, bAdExSN
        v_rAdEx = v_rAdExSNTypeII
        I0Folder = I0FolderAdExTypeIISN
        if ThetaInput:
            CSVFolder = CSVFolderAdExTypeIISNTheta
            FiguresFolder = FiguresFolderAdExTypeIISNTheta
        elif VariableDrive:
            CSVFolder = CSVFolderAdExTypeIISNVarDrive
            FiguresFolder = FiguresFolderAdExTypeIISNVarDrive	
        else:
            CSVFolder = CSVFolderAdExTypeIISN
            FiguresFolder = FiguresFolderAdExTypeIISN
    # AdEx Type II with Andronov Hopf bifurcation
    elif ChosenSubmodelAdEx==2:
        glAdEx, CAdEx, ElAdEx, DeltaTAdEx, V_TAdEx, tau_wAdEx, aAdEx, bAdEx = glAdExAH, CAdExAH, ElAdExAH, DeltaTAdExAH, V_TAdExAH, tau_wAdExAH, aAdExAH, bAdExAH
        v_rAdEx = v_rAdExAH
        u0 = w0AdExAH # Override previous value of u0
        I0Folder = I0FolderAdExTypeIIAndrHopf	
        if ThetaInput:
            CSVFolder = CSVFolderAdExTypeIIAndrHopfTheta
            FiguresFolder = FiguresFolderAdExTypeIIAndrHopfTheta
        elif VariableDrive:
            CSVFolder = CSVFolderAdExTypeIIAndrHopfVarDrive
            FiguresFolder = FiguresFolderAdExTypeIIAndrHopfVarDrive	
        else:
            CSVFolder = CSVFolderAdExTypeIIAndrHopf 
            FiguresFolder = FiguresFolderAdExTypeIIAndrHopf
    Cm = CAdEx
    v_th = v_thAdEx
    v_r = v_rAdEx
# Via's model (from Via and Canavier, 2022)
elif ChosenModel==4:
    Model = GuillemModel
    threshold_cond = threshold_condGuillem
    IntegrationStep = 0.01

    I0Folder = I0FolderGuillem
    if ThetaInput:
        CSVFolder = CSVFolderGuillemTheta
        FiguresFolder = FiguresFolderGuillemTheta
    elif VariableDrive and not UniformReversal:
        if abs(Esyn/b2.mV - -75.)<1e-3:
            CSVFolder = CSVFolderGuillemVarDriveHyper
            FiguresFolder = FiguresFolderGuillemVarDriveHyper
        elif abs(Esyn/b2.mV - -55.)<1e-3:
            CSVFolder = CSVFolderGuillemVarDriveShunt
            FiguresFolder = FiguresFolderGuillemVarDriveShunt
    elif VariableDrive and UniformReversal:
        CSVFolder = CSVFolderGuillemVarDriveUniform
        FiguresFolder = FiguresFolderGuillemVarDriveUniform
    else:
        CSVFolder = CSVFolderGuillem
        FiguresFolder = FiguresFolderGuillem

# Izhikevich Type II resonator Modified
elif ChosenModel==5:
    Model = IzhiTypeI
    threshold_cond = threshold_condIzhiTypeI
    reset_cond = reset_condIzhiTypeI
    v_th = v_thIzTyI
    v_r = v_rIzTyI
    aIzTyII = aIzTyI
    uUnits = uUnitsIzTyI
    v0 = v0IzTyI
    u0 = u0IzTyI
    Cm = CmIzTyI

    I0Folder = I0FolderIzhiTypeI
    if ThetaInput:
        CSVFolder = CSVFolderIzhiTypeIIThetaMod
        FiguresFolder = FiguresFolderIzhiTypeIIThetaMod
    elif VariableDrive:
        CSVFolder = CSVFolderIzhiTypeIIVarDriveMod
        FiguresFolder = FiguresFolderIzhiTypeIIVarDriveMod
    else:
        CSVFolder = CSVFolderIzhiTypeIIMod
        FiguresFolder = FiguresFolderIzhiTypeIIMod

# Wang Buzsaki's model (from Wang Buzsaki 1996)
elif ChosenModel==6:
    Model = WangBuzsakiModel
    threshold_cond = threshold_condWangBuzsaki
    IntegrationStep = 0.01

    I0Folder = I0FolderGuillem
    if ThetaInput:
        CSVFolder = "."#CSVFolderGuillemTheta
        FiguresFolder = "."#FiguresFolderGuillemTheta
    elif VariableDrive:
        if abs(Esyn/b2.mV - -75.)<1e-3:
            CSVFolder = CSVFolderWBVarDriveHyper
            FiguresFolder = FiguresFolderWBVarDriveHyper
        else:
            CSVFolder = CSVFolderWBVarDriveShunt
            FiguresFolder = FiguresFolderWBVarDriveShunt	
    else:
        CSVFolder = "."#CSVFolderGuillem
        FiguresFolder = "."#FiguresFolderGuillem
