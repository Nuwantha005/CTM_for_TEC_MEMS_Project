' ============================================================================
' RadialTEC_5Stage.swp - 5-Stage Radial Thermoelectric Cooler Model
' ============================================================================
' This macro creates a fully parameterized 5-stage radial TEC for MEMS cooling
' Based on the thermal network model for chip cooling applications
'
' Author: Generated for MEMS Project
' Date: November 2025
'
' Key Features:
' - 5 TEC stages with radial expansion factor (k_r)
' - P-type and N-type thermoelectric legs
' - Copper interconnects and outerconnects
' - Inter-stage ceramic insulators
' - Azimuthal insulation gaps
' - Central heat spreader cylinder
'
' Units: All dimensions in micrometers (um) internally, converted to mm for SW
' ============================================================================

Option Explicit

' ============================================================================
' GLOBAL PARAMETERS - DESIGN VARIABLES
' ============================================================================
' Geometric Parameters
Const N_STAGES As Integer = 5                   ' Number of TEC stages
Const THETA_DEG As Double = 30                  ' Wedge angle (degrees)
Const R_CYL As Double = 1000                    ' Inner cylinder radius (um)
Const W_CHIP As Double = 10000                  ' Chip width/length (um)

' Radial Expansion
Const K_R As Double = 1.3                       ' Radial expansion factor L(i+1)/L(i)

' Insulation Widths
Const W_IS As Double = 40                       ' Inter-stage insulation width (um)
Const W_AZ As Double = 20                       ' Azimuthal insulation width (um)

' Layer Thicknesses
Const T_TE As Double = 50                       ' TE leg thickness (um)
Const T_CHIP As Double = 50                     ' Chip layer thickness (um)
Const T_SOI As Double = 100                     ' SOI layer thickness (um)

' Interconnect Parameters (as ratios of stage length)
Const ALPHA_IC As Double = 0.15                 ' Interconnect ratio w_ic/L_i
Const ALPHA_OC As Double = 0.15                 ' Outerconnect ratio w_oc/L_i
Const BETA_IC_DEG As Double = 5                 ' Interconnect azimuthal angle (deg)
Const BETA_OC_DEG As Double = 5                 ' Outerconnect azimuthal angle (deg)
Const T_IC As Double = 20                       ' Interconnect thickness (um)
Const T_OC As Double = 20                       ' Outerconnect thickness (um)

' TSV Parameters
Const R_TSV As Double = 10                      ' TSV radius (um)
Const P_TSV As Double = 20                      ' TSV pitch (um)

' ============================================================================
' MATERIAL CONSTANTS
' ============================================================================
' Colors for visualization (RGB values)
Const COLOR_P_TYPE As Long = &HFF0000           ' Blue for P-type
Const COLOR_N_TYPE As Long = &H0000FF           ' Red for N-type
Const COLOR_COPPER As Long = &H00BFFF           ' Copper color
Const COLOR_INSULATOR As Long = &HC0C0C0        ' Gray for insulators
Const COLOR_SILICON As Long = &H808080          ' Dark gray for silicon

' ============================================================================
' CALCULATED PARAMETERS (Derived from design variables)
' ============================================================================
Private Pi As Double
Private ThetaRad As Double
Private R_Base As Double
Private L_TotalActive As Double
Private GeoSum As Double
Private L(1 To 5) As Double                     ' Stage lengths
Private R_Start(1 To 5) As Double               ' Stage start radii
Private R_End(1 To 5) As Double                 ' Stage end radii
Private W_IC(1 To 5) As Double                  ' Interconnect widths per stage
Private W_OC(1 To 5) As Double                  ' Outerconnect widths per stage

' SolidWorks Objects
Private swApp As SldWorks.SldWorks
Private swModel As SldWorks.ModelDoc2
Private swPart As SldWorks.PartDoc
Private swSketchMgr As SldWorks.SketchManager
Private swFeatureMgr As SldWorks.FeatureManager
Private swSelMgr As SldWorks.SelectionMgr

' ============================================================================
' MAIN ENTRY POINT
' ============================================================================
Sub Main()
    ' Initialize SolidWorks
    Set swApp = Application.SldWorks
    
    ' Calculate all derived parameters
    Call CalculateParameters
    
    ' Create new part document
    Call CreateNewPart
    
    ' Add equations to the model for parametric control
    Call AddEquationsToModel
    
    ' Create the geometry layer by layer
    Call CreateCentralCylinder
    Call CreateAllTECStages
    Call CreateInterStageInsulators
    
    ' Rebuild and zoom to fit
    swModel.EditRebuild3
    swModel.ViewZoomtofit2
    
    MsgBox "5-Stage Radial TEC Model Created Successfully!" & vbCrLf & _
           "Total Stages: " & N_STAGES & vbCrLf & _
           "Wedge Angle: " & THETA_DEG & " degrees" & vbCrLf & _
           "Expansion Factor: " & K_R, vbInformation, "RadialTEC Generator"
End Sub

' ============================================================================
' PARAMETER CALCULATION
' ============================================================================
Private Sub CalculateParameters()
    Dim i As Integer
    
    ' Constants
    Pi = 3.14159265358979
    ThetaRad = THETA_DEG * Pi / 180
    
    ' Base radius (covers square chip)
    R_Base = W_CHIP / Sqr(2)
    
    ' Available length for TE material
    ' L_total_active = (R_base - R_cyl) - (N_stages + 1) * w_is
    L_TotalActive = (R_Base - R_CYL) - (N_STAGES + 1) * W_IS
    
    ' Geometric series sum for k_r
    ' geo_sum = (1 - k_r^N) / (1 - k_r)
    If Abs(K_R - 1) < 0.0001 Then
        GeoSum = N_STAGES  ' Handle k_r = 1 case
    Else
        GeoSum = (1 - K_R ^ N_STAGES) / (1 - K_R)
    End If
    
    ' Calculate stage lengths: L_i = L_1 * k_r^(i-1)
    L(1) = L_TotalActive / GeoSum
    For i = 2 To N_STAGES
        L(i) = L(1) * K_R ^ (i - 1)
    Next i
    
    ' Calculate stage radii
    R_Start(1) = R_CYL + W_IS
    R_End(1) = R_Start(1) + L(1)
    For i = 2 To N_STAGES
        R_Start(i) = R_End(i - 1) + W_IS
        R_End(i) = R_Start(i) + L(i)
    Next i
    
    ' Calculate interconnect/outerconnect widths per stage
    For i = 1 To N_STAGES
        W_IC(i) = ALPHA_IC * L(i)
        W_OC(i) = ALPHA_OC * L(i)
    Next i
    
    ' Debug output
    Debug.Print "=== 5-Stage Radial TEC Parameters ==="
    Debug.Print "R_Base: " & R_Base & " um"
    Debug.Print "L_TotalActive: " & L_TotalActive & " um"
    Debug.Print "GeoSum: " & GeoSum
    For i = 1 To N_STAGES
        Debug.Print "Stage " & i & ": L=" & Format(L(i), "0.00") & " um, " & _
                    "r_start=" & Format(R_Start(i), "0.00") & " um, " & _
                    "r_end=" & Format(R_End(i), "0.00") & " um"
    Next i
End Sub

' ============================================================================
' CREATE NEW PART DOCUMENT
' ============================================================================
Private Sub CreateNewPart()
    Dim templatePath As String
    
    ' Get default part template
    templatePath = swApp.GetUserPreferenceStringValue(swUserPreferenceStringValue_e.swDefaultTemplatePart)
    
    If templatePath = "" Then
        templatePath = "C:\ProgramData\SolidWorks\SOLIDWORKS 2024\templates\Part.prtdot"
    End If
    
    ' Create new part
    Set swModel = swApp.NewDocument(templatePath, 0, 0, 0)
    Set swPart = swModel
    Set swSketchMgr = swModel.SketchManager
    Set swFeatureMgr = swModel.FeatureManager
    Set swSelMgr = swModel.SelectionManager
    
    ' Set document units to micrometers (we'll work in mm and scale)
    ' SolidWorks uses meters internally, so we convert um to m
    ' 1 um = 0.001 mm = 0.000001 m
End Sub

' ============================================================================
' ADD EQUATIONS TO MODEL
' ============================================================================
Private Sub AddEquationsToModel()
    Dim swEquationMgr As SldWorks.EquationMgr
    Dim eqIdx As Integer
    Dim i As Integer
    
    Set swEquationMgr = swModel.GetEquationMgr
    
    ' Add global variables as equations
    eqIdx = swEquationMgr.Add2(-1, """N_stages"" = " & N_STAGES, True)
    eqIdx = swEquationMgr.Add2(-1, """theta_deg"" = " & THETA_DEG, True)
    eqIdx = swEquationMgr.Add2(-1, """R_cyl"" = " & R_CYL & "um", True)
    eqIdx = swEquationMgr.Add2(-1, """w_chip"" = " & W_CHIP & "um", True)
    eqIdx = swEquationMgr.Add2(-1, """k_r"" = " & K_R, True)
    eqIdx = swEquationMgr.Add2(-1, """w_is"" = " & W_IS & "um", True)
    eqIdx = swEquationMgr.Add2(-1, """w_az"" = " & W_AZ & "um", True)
    eqIdx = swEquationMgr.Add2(-1, """t_TE"" = " & T_TE & "um", True)
    eqIdx = swEquationMgr.Add2(-1, """alpha_ic"" = " & ALPHA_IC, True)
    eqIdx = swEquationMgr.Add2(-1, """alpha_oc"" = " & ALPHA_OC, True)
    
    ' Add calculated values
    For i = 1 To N_STAGES
        eqIdx = swEquationMgr.Add2(-1, """L_" & i & """ = " & Format(L(i), "0.00") & "um", True)
        eqIdx = swEquationMgr.Add2(-1, """r_start_" & i & """ = " & Format(R_Start(i), "0.00") & "um", True)
        eqIdx = swEquationMgr.Add2(-1, """r_end_" & i & """ = " & Format(R_End(i), "0.00") & "um", True)
    Next i
End Sub

' ============================================================================
' CREATE CENTRAL HEAT SPREADER CYLINDER
' ============================================================================
Private Sub CreateCentralCylinder()
    Dim swSketch As SldWorks.Sketch
    Dim swFeat As SldWorks.Feature
    Dim radiusMM As Double
    Dim heightMM As Double
    
    ' Convert to mm (SolidWorks default unit)
    radiusMM = R_CYL / 1000
    heightMM = T_TE / 1000
    
    ' Create sketch on Top plane
    swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    swSketchMgr.InsertSketch True
    
    ' Draw circle for central cylinder
    swSketchMgr.CreateCircle 0, 0, 0, radiusMM, 0, 0
    
    ' Exit sketch and extrude
    swSketchMgr.InsertSketch True
    
    ' Boss extrude
    swModel.Extension.SelectByID2 "Sketch1", "SKETCH", 0, 0, 0, False, 0, Nothing, 0
    swFeatureMgr.FeatureExtrusion2 True, False, False, 0, 0, heightMM, 0, _
        False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
    
    ' Rename feature
    swModel.Extension.SelectByID2 "Boss-Extrude1", "BODYFEATURE", 0, 0, 0, False, 0, Nothing, 0
    swModel.SelectedFeatureProperties 0, 0, 0, 0, 0, 0, 0, True, False, "CentralCylinder"
End Sub

' ============================================================================
' CREATE ALL TEC STAGES
' ============================================================================
Private Sub CreateAllTECStages()
    Dim stage As Integer
    
    For stage = 1 To N_STAGES
        Call CreateTECStage(stage)
        Debug.Print "Created Stage " & stage
    Next stage
End Sub

' ============================================================================
' CREATE SINGLE TEC STAGE
' ============================================================================
Private Sub CreateTECStage(stageNum As Integer)
    Dim swSketch As SldWorks.Sketch
    Dim swFeat As SldWorks.Feature
    
    ' Create P-type leg
    Call CreateTELeg(stageNum, "P")
    
    ' Create N-type leg
    Call CreateTELeg(stageNum, "N")
    
    ' Create copper interconnect at inner radius
    Call CreateInterconnect(stageNum)
    
    ' Create copper outerconnect at outer radius
    Call CreateOuterconnect(stageNum)
End Sub

' ============================================================================
' CREATE THERMOELECTRIC LEG (P or N type)
' ============================================================================
Private Sub CreateTELeg(stageNum As Integer, legType As String)
    Dim r1 As Double, r2 As Double
    Dim theta1 As Double, theta2 As Double
    Dim heightMM As Double
    Dim sketchName As String
    Dim featName As String
    Dim halfAngle As Double
    Dim azGapRad As Double
    
    ' Convert dimensions to mm
    r1 = (R_Start(stageNum) + W_IC(stageNum)) / 1000    ' After interconnect
    r2 = (R_End(stageNum) - W_OC(stageNum)) / 1000       ' Before outerconnect
    heightMM = T_TE / 1000
    
    ' Calculate azimuthal gap in radians
    azGapRad = (W_AZ / 1000) / ((r1 + r2) / 2)  ' Arc length / avg radius
    halfAngle = ThetaRad / 4  ' Half of half-wedge for P or N leg
    
    ' Determine angular position based on leg type
    If legType = "P" Then
        ' P-type on one side of centerline
        theta1 = azGapRad / 2
        theta2 = ThetaRad / 2 - azGapRad / 2
        featName = "Stage" & stageNum & "_P_Leg"
    Else
        ' N-type on other side of centerline (negative angles)
        theta1 = -ThetaRad / 2 + azGapRad / 2
        theta2 = -azGapRad / 2
        featName = "Stage" & stageNum & "_N_Leg"
    End If
    
    ' Create sketch on Top plane
    swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    swSketchMgr.InsertSketch True
    
    ' Draw wedge-shaped leg
    Call DrawWedgeShape(r1, r2, theta1, theta2)
    
    ' Exit sketch and extrude
    swSketchMgr.InsertSketch True
    
    ' Get the latest sketch
    Dim sketchCount As Integer
    sketchCount = swModel.GetFeatureCount
    
    ' Boss extrude
    swFeatureMgr.FeatureExtrusion2 True, False, False, 0, 0, heightMM, 0, _
        False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
    
    ' Rename feature
    Call RenameLastFeature(featName)
End Sub

' ============================================================================
' CREATE COPPER INTERCONNECT (Inner connection between P and N)
' ============================================================================
Private Sub CreateInterconnect(stageNum As Integer)
    Dim r1 As Double, r2 As Double
    Dim theta1 As Double, theta2 As Double
    Dim heightMM As Double
    Dim betaRad As Double
    Dim featName As String
    
    ' Convert dimensions to mm
    r1 = R_Start(stageNum) / 1000
    r2 = (R_Start(stageNum) + W_IC(stageNum)) / 1000
    heightMM = T_IC / 1000
    
    ' Interconnect spans the azimuthal angle
    betaRad = BETA_IC_DEG * Pi / 180
    theta1 = -betaRad / 2
    theta2 = betaRad / 2
    
    featName = "Stage" & stageNum & "_Interconnect"
    
    ' Create sketch on Top plane
    swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    swSketchMgr.InsertSketch True
    
    ' Draw interconnect wedge
    Call DrawWedgeShape(r1, r2, theta1, theta2)
    
    ' Exit sketch and extrude
    swSketchMgr.InsertSketch True
    
    ' Boss extrude (offset from TE layer base)
    Dim offsetMM As Double
    offsetMM = (T_TE - T_IC) / 2 / 1000  ' Center the interconnect in TE thickness
    
    swFeatureMgr.FeatureExtrusion2 True, False, False, 0, 0, heightMM, 0, _
        False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
    
    ' Rename feature
    Call RenameLastFeature(featName)
End Sub

' ============================================================================
' CREATE COPPER OUTERCONNECT (Outer connection to next stage)
' ============================================================================
Private Sub CreateOuterconnect(stageNum As Integer)
    Dim r1 As Double, r2 As Double
    Dim theta1 As Double, theta2 As Double
    Dim heightMM As Double
    Dim betaRad As Double
    Dim featName As String
    
    ' Convert dimensions to mm
    r1 = (R_End(stageNum) - W_OC(stageNum)) / 1000
    r2 = R_End(stageNum) / 1000
    heightMM = T_OC / 1000
    
    ' Outerconnect spans the azimuthal angle
    betaRad = BETA_OC_DEG * Pi / 180
    theta1 = -betaRad / 2
    theta2 = betaRad / 2
    
    featName = "Stage" & stageNum & "_Outerconnect"
    
    ' Create sketch on Top plane
    swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    swSketchMgr.InsertSketch True
    
    ' Draw outerconnect wedge
    Call DrawWedgeShape(r1, r2, theta1, theta2)
    
    ' Exit sketch and extrude
    swSketchMgr.InsertSketch True
    
    ' Boss extrude
    swFeatureMgr.FeatureExtrusion2 True, False, False, 0, 0, heightMM, 0, _
        False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
    
    ' Rename feature
    Call RenameLastFeature(featName)
End Sub

' ============================================================================
' CREATE INTER-STAGE INSULATORS
' ============================================================================
Private Sub CreateInterStageInsulators()
    Dim i As Integer
    
    ' Create insulator between central cylinder and stage 1
    Call CreateInsulatorRing(0)
    
    ' Create insulators between stages
    For i = 1 To N_STAGES - 1
        Call CreateInsulatorRing(i)
    Next i
    
    ' Create outer insulator after stage N
    Call CreateInsulatorRing(N_STAGES)
End Sub

' ============================================================================
' CREATE SINGLE INSULATOR RING
' ============================================================================
Private Sub CreateInsulatorRing(afterStage As Integer)
    Dim r1 As Double, r2 As Double
    Dim heightMM As Double
    Dim featName As String
    
    heightMM = T_TE / 1000
    
    If afterStage = 0 Then
        ' Insulator between center and stage 1
        r1 = R_CYL / 1000
        r2 = (R_CYL + W_IS) / 1000
        featName = "Insulator_Center"
    ElseIf afterStage = N_STAGES Then
        ' Outer insulator
        r1 = R_End(N_STAGES) / 1000
        r2 = (R_End(N_STAGES) + W_IS) / 1000
        featName = "Insulator_Outer"
    Else
        ' Insulator between stages
        r1 = R_End(afterStage) / 1000
        r2 = (R_End(afterStage) + W_IS) / 1000
        featName = "Insulator_S" & afterStage & "_S" & (afterStage + 1)
    End If
    
    ' Create sketch on Top plane
    swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    swSketchMgr.InsertSketch True
    
    ' Draw insulator ring (full wedge angle)
    Call DrawWedgeShape(r1, r2, -ThetaRad / 2, ThetaRad / 2)
    
    ' Exit sketch and extrude
    swSketchMgr.InsertSketch True
    
    ' Boss extrude
    swFeatureMgr.FeatureExtrusion2 True, False, False, 0, 0, heightMM, 0, _
        False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
    
    ' Rename feature
    Call RenameLastFeature(featName)
End Sub

' ============================================================================
' DRAW WEDGE SHAPE (Arc segment)
' ============================================================================
Private Sub DrawWedgeShape(r1 As Double, r2 As Double, theta1 As Double, theta2 As Double)
    Dim x1 As Double, y1 As Double
    Dim x2 As Double, y2 As Double
    Dim x3 As Double, y3 As Double
    Dim x4 As Double, y4 As Double
    
    ' Calculate corner points
    ' Point 1: Inner radius, start angle
    x1 = r1 * Cos(theta1)
    y1 = r1 * Sin(theta1)
    
    ' Point 2: Outer radius, start angle
    x2 = r2 * Cos(theta1)
    y2 = r2 * Sin(theta1)
    
    ' Point 3: Outer radius, end angle
    x3 = r2 * Cos(theta2)
    y3 = r2 * Sin(theta2)
    
    ' Point 4: Inner radius, end angle
    x4 = r1 * Cos(theta2)
    y4 = r1 * Sin(theta2)
    
    ' Draw the wedge shape
    ' Line 1: Inner radial edge (point 1 to point 2)
    swSketchMgr.CreateLine x1, y1, 0, x2, y2, 0
    
    ' Arc 1: Outer arc (point 2 to point 3)
    swSketchMgr.CreateArc 0, 0, 0, x2, y2, 0, x3, y3, 0, 1
    
    ' Line 2: Outer radial edge (point 3 to point 4)
    swSketchMgr.CreateLine x3, y3, 0, x4, y4, 0
    
    ' Arc 2: Inner arc (point 4 to point 1)
    swSketchMgr.CreateArc 0, 0, 0, x4, y4, 0, x1, y1, 0, -1
End Sub

' ============================================================================
' RENAME LAST CREATED FEATURE
' ============================================================================
Private Sub RenameLastFeature(newName As String)
    Dim swFeat As SldWorks.Feature
    Dim featCount As Integer
    
    ' Get last feature
    Set swFeat = swModel.FirstFeature
    Do While Not swFeat.GetNextFeature Is Nothing
        Set swFeat = swFeat.GetNextFeature
    Loop
    
    If Not swFeat Is Nothing Then
        swFeat.Name = newName
    End If
End Sub

' ============================================================================
' UTILITY: Convert um to mm
' ============================================================================
Private Function UmToMm(valueUm As Double) As Double
    UmToMm = valueUm / 1000
End Function

' ============================================================================
' UTILITY: Convert degrees to radians
' ============================================================================
Private Function DegToRad(valueDeg As Double) As Double
    DegToRad = valueDeg * 3.14159265358979 / 180
End Function

' ============================================================================
' CREATE TSV ARRAY (Optional - for thermal vias in evaporator zone)
' ============================================================================
Public Sub CreateTSVArray(stageNum As Integer)
    Dim r_ic As Double
    Dim numTSV As Integer
    Dim pitch As Double
    Dim i As Integer
    Dim angle As Double
    Dim x As Double, y As Double
    Dim heightMM As Double
    Dim radiusMM As Double
    
    ' Only create TSVs in evaporator zone (first 2-3 stages typically)
    If stageNum > 3 Then Exit Sub
    
    ' Position TSVs at interconnect radius
    r_ic = (R_Start(stageNum) + W_IC(stageNum) / 2) / 1000
    radiusMM = R_TSV / 1000
    heightMM = T_SOI / 1000
    pitch = P_TSV / 1000
    
    ' Calculate number of TSVs based on arc length and pitch
    Dim arcLength As Double
    arcLength = r_ic * BETA_IC_DEG * Pi / 180
    numTSV = Int(arcLength / pitch)
    
    If numTSV < 1 Then numTSV = 1
    
    ' Create TSV circles
    For i = 1 To numTSV
        angle = -BETA_IC_DEG * Pi / 360 + (i - 0.5) * (BETA_IC_DEG * Pi / 180) / numTSV
        x = r_ic * Cos(angle)
        y = r_ic * Sin(angle)
        
        ' Create sketch for each TSV
        swModel.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
        swSketchMgr.InsertSketch True
        swSketchMgr.CreateCircle x, y, 0, x + radiusMM, y, 0
        swSketchMgr.InsertSketch True
        
        ' Extrude TSV (downward into chip)
        swFeatureMgr.FeatureExtrusion2 True, False, False, 1, 0, heightMM, 0, _
            False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False
        
        Call RenameLastFeature("TSV_S" & stageNum & "_" & i)
    Next i
End Sub

' ============================================================================
' GENERATE PARAMETER SUMMARY
' ============================================================================
Public Sub PrintParameterSummary()
    Dim i As Integer
    
    Debug.Print ""
    Debug.Print "=============================================="
    Debug.Print "5-STAGE RADIAL TEC PARAMETER SUMMARY"
    Debug.Print "=============================================="
    Debug.Print ""
    Debug.Print "=== GLOBAL PARAMETERS ==="
    Debug.Print "Number of stages (N_stages): " & N_STAGES
    Debug.Print "Wedge angle (theta): " & THETA_DEG & " deg"
    Debug.Print "Inner cylinder radius (R_cyl): " & R_CYL & " um"
    Debug.Print "Chip dimension (w_chip): " & W_CHIP & " um"
    Debug.Print "Base radius (R_base): " & Format(R_Base, "0.00") & " um"
    Debug.Print ""
    Debug.Print "=== EXPANSION PARAMETERS ==="
    Debug.Print "Radial expansion factor (k_r): " & K_R
    Debug.Print "Total active length: " & Format(L_TotalActive, "0.00") & " um"
    Debug.Print ""
    Debug.Print "=== INSULATION PARAMETERS ==="
    Debug.Print "Inter-stage width (w_is): " & W_IS & " um"
    Debug.Print "Azimuthal width (w_az): " & W_AZ & " um"
    Debug.Print ""
    Debug.Print "=== THICKNESS PARAMETERS ==="
    Debug.Print "TE thickness (t): " & T_TE & " um"
    Debug.Print "Chip thickness (t_chip): " & T_CHIP & " um"
    Debug.Print "SOI thickness (t_SOI): " & T_SOI & " um"
    Debug.Print "Interconnect thickness (t_ic): " & T_IC & " um"
    Debug.Print "Outerconnect thickness (t_oc): " & T_OC & " um"
    Debug.Print ""
    Debug.Print "=== INTERCONNECT RATIOS ==="
    Debug.Print "Interconnect ratio (alpha_ic): " & ALPHA_IC
    Debug.Print "Outerconnect ratio (alpha_oc): " & ALPHA_OC
    Debug.Print "Interconnect angle (beta_ic): " & BETA_IC_DEG & " deg"
    Debug.Print "Outerconnect angle (beta_oc): " & BETA_OC_DEG & " deg"
    Debug.Print ""
    Debug.Print "=== STAGE-BY-STAGE GEOMETRY ==="
    Debug.Print "Stage | Length (um) | r_start (um) | r_end (um) | w_ic (um) | w_oc (um)"
    Debug.Print "----------------------------------------------------------------------"
    For i = 1 To N_STAGES
        Debug.Print Format(i, "  0  ") & " | " & _
                    Format(L(i), "0000.00") & "    | " & _
                    Format(R_Start(i), "0000.00") & "     | " & _
                    Format(R_End(i), "0000.00") & "   | " & _
                    Format(W_IC(i), "000.00") & "   | " & _
                    Format(W_OC(i), "000.00")
    Next i
    Debug.Print "=============================================="
End Sub
