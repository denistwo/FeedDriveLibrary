within ;
package FeedDriveLibrary
  import SI = Modelica.SIunits;
  import Modelica.Constants;

  package Information
    extends Modelica.Icons.Information;
     annotation (Documentation(info="<html>
<p>
The objective of the Feed Drive Library to model linear axes in machine tools and production machines. Therefore the library contains key elements such as converter, motor, gear, clutch etc. 
The models are designed in a way so that they can easily be parameterized
</p>

The objective of the Feed Drive Library is to model linear axes in machine tools and production machines. Therefore, the library contains key elements such as converter, motor, gear, clutch etc. 
The models are designed in a way so that they can easily be parameterized with typical vendor datasheets. Thus, we defined a simpler motor model in comparison to the models in electric machine models in the Modelica Standard Library. Further, the models do not only contain the behavior equations, but also comprise the metrics and requirements that are important to choose adequate components for designing feed drive axes. These requirements can also be used for system optimization.
The library is based on basic elements from the Standard Library and from the Optimization Library from Dymola. Regarding the Standard Library these are icons and different adapted models. Regarding the Optimization Library we had a look at the Criteria models (mean, max, etc) and changed these models for our specific purposes.

<br \><br \>
<h4>
Disclaimer
</h4>
No liability can be accepted for any errors or omissions.
<h4>
Reference
</h4>
When using the library please cite 
Özdemir, D.; Motschke, T.; Herfs, W.; Brecher, C.: Modelica Library for Feed Drive Systems. In: 11th International Modelica Conference, Paris, 2015, pp. 117-125
<h4>
License
</h4>
The library is based on the Modelica Standard Library. For our own work MSD/MIT shall apply. 
<h4>
Contact
</h4>
Laboratory for Machine Tools and Production Systems<br \>
RWTH Aachen University<br \>
<a href='http://www.wzl.rwth-aachen.de' >Website of WZL Aachen</a>
<h4>
Acknowledgement
</h4>
We gratefully acknowledge funding from the German Research Foundation (DFG) in the project &ldquo;Optimierung des Systementwurfs von Maschinen und Anlagen auf Basis komponentenorientierter Verhaltensmodelle&rdquo;
</html>
"));

  end Information;

  package Examples

    model S1_Operation_Siemens

      Modelica.Mechanics.Rotational.Sources.Torque torque annotation(Placement(transformation(extent={{40,-50},
                {60,-30}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T = 313.15) annotation(Placement(transformation(extent={{-60,-50},
                {-40,-30}})));
      FeedDriveMotor.Motors.PMSM Siemens_1FT7086_AC7(n_n = 2000, n_polePairs = 5, M_n = 22.5, I_n = 9.2, M_0_100 = 28, I_0_100 = 10.6, M_0_60 = 23, I_0_60 = 8.6, J = 0.00636, n_max = 8000, M_max = 120, I_max = 54, R_ph_20 = 0.46, L_D = 0.0085,
        T_max=413.15,
        U_eff_max=380,
        T_th=3600,
        T_start=313.15)                                                                                                     annotation(Placement(transformation(extent={{-62,-12},
                {-8,14}})));
      Modelica.Blocks.Sources.Constant const(k=93)
        annotation (Placement(transformation(extent={{-92,-30},{-72,-10}})));
      Modelica.Blocks.Sources.Constant const1(k=10)
        annotation (Placement(transformation(extent={{10,-50},{30,-30}})));
    equation
      connect(fixedTemperature.port, Siemens_1FT7086_AC7.heatPort) annotation (Line(
          points={{-40,-40},{-21.5,-40},{-21.5,-10.7}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(torque.flange, Siemens_1FT7086_AC7.flange) annotation (Line(
          points={{60,-40},{80,-40},{80,1},{-9.35,1}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(torque.tau, const1.y) annotation (Line(
          points={{38,-40},{31,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(const.y, Siemens_1FT7086_AC7.u_q) annotation (Line(
          points={{-71,-20},{-48.5,-20},{-48.5,-12.26}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                -100},{100,100}}),                                                                           graphics));
    end S1_Operation_Siemens;

    model SimpleApplicationScenario_LM_dyn
      FeedDriveMotor.Motors.LM_cooling S_1FN3450_3WE(F_N = 2895, I_n = 50.7, m = 22.6, F_max = 7760, I_max = 149.6, U_eff_max = 425, k_F_0_20 = 57, R_ph_20 = 0.2, L_D = 0.0024, tau_p = 0.023, T_th = 120, T_start = 393.15, Q_N = 7.5e-05, Delta_p_N = 65000) annotation(Placement(transformation(extent = {{-50, -2}, {-10, 18}})));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 47.4) annotation(Placement(transformation(extent = {{8, -2}, {28, 18}})));
      TranslationalComponents.Force_PT1 force_PT1_1(T = 0.00001) annotation(Placement(transformation(extent = {{58, -2}, {38, 18}})));
      Cooling.DisplacementPump_Constant displacementPump(Q = 7.5e-05, T = 293.15) annotation(Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = 270, origin = {-20, 48})));
      Cooling.Tank        tank_Simple annotation(Placement(transformation(extent = {{-10, 32}, {10, 52}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor meas annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {8, 28})));
      Modelica.Blocks.Continuous.PI PI_I(T = 0.012, k = 50) annotation(Placement(transformation(extent = {{-92, -22}, {-72, -2}})));
      HelpBlocks.Feedback_mirror I_feedback annotation(Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 180, origin = {-110, -12})));
      Modelica.Blocks.Sources.Constant const1(k = 0) annotation(Placement(transformation(extent = {{-180, -22}, {-160, -2}})));
      Modelica.Blocks.Sources.Step step(startTime = 0.005, height = -100) annotation(Placement(transformation(extent = {{38, 26}, {58, 46}})));
      Modelica.Blocks.Nonlinear.FixedDelay fixedDelay1(delayTime = 62.5e-6) annotation(Placement(transformation(extent = {{-70, 18}, {-90, 38}})));
      FeedDriveMotor.Components.v_p_controller v_p_controller(T_nn = 0.003, delay = 62.5e-6, K_V = 40, K_pn = 1400) annotation(Placement(transformation(extent = {{-146, -22}, {-126, -2}})));
    equation
      connect(force_PT1_1.flange, mass.flange_b) annotation(Line(points = {{38, 8}, {28, 8}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(S_1FN3450_3WE.coolingPort_a, displacementPump.port_b) annotation(Line(points = {{-25.8, 13.8}, {-30, 13.8}, {-30, 48}, {-28, 48}}, color = {0, 0, 0}, pattern = LinePattern.Solid, thickness = 0.25, smooth = Smooth.None));
      connect(displacementPump.port_a, tank_Simple.port_a) annotation(Line(points = {{-12, 48}, {0, 48}}, color = {0, 0, 0}, pattern = LinePattern.Solid, thickness = 0.25, smooth = Smooth.None));
      connect(mass.flange_a, S_1FN3450_3WE.f) annotation(Line(points = {{8, 8}, {-11, 8}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(S_1FN3450_3WE.coolingPort_b, tank_Simple.port_a) annotation(Line(points = {{-14.2, 13.8}, {-10, 13.8}, {-10, 48}, {0, 48}}, color = {0, 0, 0}, pattern = LinePattern.Solid, smooth = Smooth.None));
      connect(meas.flange, mass.flange_a) annotation(Line(points = {{8, 18}, {8, 8}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(PI_I.y, S_1FN3450_3WE.u_q) annotation(Line(points = {{-71, -12}, {-40, -12}, {-40, -2.2}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(I_feedback.y, PI_I.u) annotation(Line(points = {{-101, -12}, {-94, -12}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(S_1FN3450_3WE.i_q, fixedDelay1.u) annotation(Line(points = {{-40, 17}, {-40, 28}, {-68, 28}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(fixedDelay1.y, I_feedback.u2) annotation(Line(points = {{-91, 28}, {-110, 28}, {-110, -4}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(const1.y, v_p_controller.x_s) annotation(Line(points = {{-159, -12}, {-148, -12}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(v_p_controller.i_s, I_feedback.u1) annotation(Line(points = {{-125, -12}, {-118, -12}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(meas.v, v_p_controller.v_a) annotation(Line(points = {{8, 39}, {8, 58}, {-136, 58}, {-136, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(step.y, force_PT1_1.f) annotation(Line(points = {{59, 36}, {72, 36}, {72, 8}, {60, 8}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}})));
    end SimpleApplicationScenario_LM_dyn;

    model h_i_opt

      FeedDriveMotor.Motors.PMSM PSM(n_n = 2000, n_polePairs = 5, M_n = 22.5, I_n = 9.2, M_0_100 = 28, I_0_100 = 10.6, M_0_60 = 23, I_0_60 = 8.6, J = 0.00636, n_max = 8000, M_max = 120, I_max = 54, U_eff_max = 425, L_D = 0.0085, C_F = 0, R_ph_20 = 0.001, T_max = 413.15, T_th = 3600, T_start = 393.15) annotation(Placement(transformation(extent = {{-84, -10}, {-44, 10}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature T_U(T = 313.15) annotation(Placement(transformation(extent={{-80,-60},
                {-60,-40}})));
      Modelica.Mechanics.Translational.Components.Mass table(m = 2500) annotation(Placement(transformation(extent = {{30, -10}, {50, 10}})));
      Modelica.Blocks.Sources.Ramp ramp(height=55, duration=0.3)      annotation(Placement(transformation(extent={{-116,
                -30},{-96,-10}})));
      LinearActuators.BallScrew.BallScrewDrive BSD(mu = 0, J = 0.003946 / 2, rho = 0, bearingConfiguration = LinearActuators.BallScrew.BearingConfiguration.fixed_free, l_S = 1.4, C_am = 1.5e5, d = 0.063, P = 0.03) annotation(Placement(transformation(extent = {{4, -10}, {24, 10}})));
      Modelica.Mechanics.Translational.Components.Fixed machineBed annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-10, -30})));
      RotationalComponents.TimingBelt.PulleyConstraints TBD_small(d_eff = 0.08149, z = 32, J = 0.000655, eta = 1, exponent = 0.827963, factor = 19.47) annotation(Placement(transformation(extent = {{-40, -10}, {-20, 10}})));
      RotationalComponents.TimingBelt.BeltPulley TBD_large(d_eff = 0.08149, z = 32, J = 0.000655) annotation(Placement(transformation(extent = {{0, -10}, {-20, 10}})));
    equation
      connect(T_U.port, PSM.heatPort) annotation(Line(points={{-60,-50},{-54,
              -50},{-54,-9}},                                                                       color = {191, 0, 0}, smooth = Smooth.None));
      connect(BSD.flangeT, table.flange_a) annotation(Line(points = {{24.2, 0}, {30, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(machineBed.flange, BSD.supportT1) annotation(Line(points = {{-10, -30}, {3.8, -30}, {3.8, -7}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(TBD_small.flangeR, PSM.flange) annotation(Line(points = {{-40, 0}, {-45, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(TBD_small.flangeT, TBD_large.flangeT) annotation(Line(points = {{-20, 0}, {-20, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(TBD_large.flangeR, BSD.flangeR) annotation(Line(points = {{0, 0}, {5, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(ramp.y, PSM.u_q) annotation (Line(
          points={{-95,-20},{-74,-20},{-74,-10.2}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                -100},{100,100}}),                                                                           graphics));
    end h_i_opt;

    model FrequencyAnalysis

      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature T_U(T = 298.15) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={82,-58})));
      Modelica.Blocks.Continuous.PI PI_w(k=8.264672852, T=0.03)
                                                             annotation(Placement(transformation(extent={{4,-68},
                {24,-48}})));
      Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = 180, origin={-4,12})));
      Modelica.Mechanics.Rotational.Components.SpringDamper encoderMounting(c = 3000, d = 0.02) annotation(Placement(transformation(extent={{44,2},{
                64,22}})));
      Modelica.Mechanics.Rotational.Components.Inertia encoder(a(start = 0), w(start = 0), J = 3e-6) annotation(Placement(transformation(extent={{14,2},{
                34,22}})));
      Modelica.Mechanics.Rotational.Components.SpringDamper screwShaft(c=2000, d=
           1)
        annotation (Placement(transformation(extent={{92,-30},{112,-10}})));
      LinearActuators.BallScrew.BallScrewDrive ballscrew(C_am = 6.20E+04, C_0am = 1.04E+05, l_S = 0.85, R_nu = 2e9, R_S = 0.6e9,                      rho = 0, mu = 0,
        d=0.063,
        P=0.02,
        J=0.01)                                                                                                     annotation(Placement(transformation(extent={{122,-30},
                {142,-10}})));
      Modelica.Mechanics.Translational.Components.Fixed machineBed annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin={112,-58})));
      HelpBlocks.to_mm to_mm annotation(Placement(transformation(extent={{18,50},
                {-2,70}})));
      Modelica.Mechanics.Translational.Sensors.PositionSensor positionSensor annotation(Placement(transformation(extent={{124,50},
                {104,70}})));
      Modelica.Mechanics.Translational.Sources.Force force
        annotation (Placement(transformation(extent={{200,10},{220,30}})));
      Modelica.Mechanics.Translational.Components.SpringDamper screwShaft1(c=1.38E8,
          d=200e3)                                                                            annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin={164,-20})));
      Modelica.Mechanics.Translational.Components.Mass machineBase(
        s(start=0),
        v(start=0),
        a(start=0),
        m=4000)     annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={208,-20})));
      Modelica.Blocks.Sources.Step step1(           startTime=0.5, height=6e3
            /(1.59/0.5))
        annotation (Placement(transformation(extent={{166,10},{186,30}})));
      FeedDriveMotor.Motors_withAutoTunedCurrentController.PMSM pMSM(U_eff_max = 425, n_n = 1500, n_polePairs = 4, M_n = 61, I_n = 20.5, M_0_100 = 70, I_0_100 = 22.3, M_0_60 = 58, I_0_60 = 18.1,             n_max = 5600, M_max = 220, I_max = 107, R_ph_20 = 0.22, L_D = 0.0052,
        J=0.260,
        T_max=413.15,
        T_th=3300,
        T_start=413.15) "1FT6108-8AB71"                                                                                     annotation(Placement(transformation(extent={{34,-30},
                {74,-10}})));
      Modelica.Blocks.Interfaces.RealInput u1 "Connector of Real input signal"
        annotation (Placement(transformation(extent={{-96,-78},{-56,-38}})));
      Modelica.Blocks.Interfaces.RealOutput w1
        "Absolute angular velocity of flange as output signal" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-60,12})));
      Modelica.Blocks.Interfaces.RealOutput y1
        "Connector of Real output signal containing input signal u in another unit"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-54,60})));
    equation
      connect(speedSensor.flange, encoder.flange_a) annotation(Line(points={{6,12},{
              14,12}},                                                                               color = {0, 0, 0}, smooth = Smooth.None));
      connect(encoder.flange_b, encoderMounting.flange_a) annotation(Line(points={{34,12},
              {44,12}},                                                                                    color = {0, 0, 0}, smooth = Smooth.None));
      connect(screwShaft.flange_b, ballscrew.flangeR) annotation (Line(
          points={{112,-20},{123,-20}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(machineBed.flange, ballscrew.supportT1) annotation (Line(
          points={{112,-58},{121.8,-58},{121.8,-27}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(to_mm.u,positionSensor. s) annotation (Line(
          points={{20,60},{103,60}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ballscrew.flangeT, screwShaft1.flange_a) annotation (Line(
          points={{142.2,-20},{154,-20}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(screwShaft1.flange_b, machineBase.flange_a) annotation (Line(
          points={{174,-20},{198,-20}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(positionSensor.flange, screwShaft1.flange_a) annotation (Line(
          points={{124,60},{148,60},{148,-20},{154,-20}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(force.flange, machineBase.flange_b) annotation (Line(
          points={{220,20},{238,20},{238,-20},{218,-20}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(force.f, step1.y) annotation (Line(
          points={{198,20},{187,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PI_w.y, pMSM.i_set) annotation (Line(
          points={{25,-58},{44,-58},{44,-30}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pMSM.heatPort, T_U.port) annotation (Line(
          points={{64,-29},{66,-29},{66,-58},{72,-58}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(pMSM.flange, screwShaft.flange_a) annotation (Line(
          points={{73,-20},{92,-20}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(encoderMounting.flange_b, screwShaft.flange_a) annotation (Line(
          points={{64,12},{82,12},{82,-20},{92,-20}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(PI_w.u, u1) annotation (Line(
          points={{2,-58},{-76,-58}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(speedSensor.w, w1) annotation (Line(
          points={{-15,12},{-60,12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(to_mm.y, y1) annotation (Line(
          points={{-3,60},{-54,60}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-250,
                -100},{250,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-250, -100}, {250, 100}})));
    end FrequencyAnalysis;
  end Examples;

  package FeedDriveMotor
    package Motors
      model PMSM
        "Model of a synchronous servo motor that can be parametrized with typical supplier data such as from Siemens 1FT6. Field weakening current is set to zero."
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_n
          "Rated speed"                                                           annotation(Dialog(group = "Engineering data"));
        parameter Integer n_polePairs "Number of pole pairs" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_n "Rated torque (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_n "Rated current (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_0_100 "Stall torque (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_0_100 "Stall current (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_0_60 "Stall torque (60 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_0_60 "Stall current (60 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Inertia J "Motor inertia" annotation(Dialog(group = "Engineering data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_max
          "Max. permissible speed (mech.)"                                                             annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_max "Maximum torque" annotation(Dialog(group = "Limiting data"));
        parameter SI.Current I_max "Maximum current" annotation(Dialog(group = "Limiting data"));
        parameter SI.Temperature T_max(displayUnit = "degC") = 413.15
          "Permissable temperature"                                                             annotation(Dialog(group = "Limiting data"));
        parameter SI.Voltage U_eff_max = 380 "Voltage limit of motor module" annotation(Dialog(group = "Limiting data"));
        parameter SI.ElectricalTorqueConstant k_T_0_100 = M_0_100 / I_0_100
          "Rated torque constant"                                                                   annotation(Dialog(group = "Physical constants"));
        parameter SI.Resistance R_ph_20 "Winding resistance at 20degC" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inductance L_D "Rotating field inductance" annotation(Dialog(group = "Physical constants"));
        parameter SI.Time T_th(displayUnit = "min") "Thermal time constant" annotation(Dialog(group = "Physical constants"));
        parameter Real C_F = (T_ref_100 / R_heat - Q_R_N) / SI.Conversions.from_rpm(n_n) ^ 1.5
          "Constant relating friction torque and square root of angular velocity"
                                                                                                              annotation(Dialog(group = "Physical constants"));
        parameter SI.Temperature T_start "Start temperature" annotation(Dialog(group = "Start conditions"));
        parameter Real S_n = 1
          "Safety factor (>=1) for rotation speed (S_n * n <= n_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_M = 1
          "Safety factor (>=1) for torque (S_n * M <= M_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_U = 1
          "Safety factor (>=1) for voltage (S_n * u <= u_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_T = 1
          "Safety factor (>=1) for temperature (S_n * T <= T_max)"                      annotation(Dialog(tab = "Safety factors"));
        Real con_n_max;
        Real con_M_max;
        Real con_T_max "End of period value, for simulation of one period";
        Real con_U_max;
        Real util_n_max;
        Real util_M_max;
        Real util_T_max;
        Real util_U_max;
        Real obj_rms_current;
        Real obj_absMax_current;
        Real obj_power_average;
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange
          "Flange of shaft"                                                        annotation(Placement(transformation(extent = {{180, -10}, {200, 10}}), iconTransformation(extent = {{180, -10}, {200, 10}})));
      protected
        parameter SI.HeatCapacity C = T_th / R_heat "Heat Capacity";
        parameter SI.HeatFlowRate Q_el_100 = 3 * I_0_100 ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at M = M_0_100";
        parameter SI.HeatFlowRate Q_R_N = 3 * (M_n / k_T_0_100) ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at nominal point";
        parameter SI.ThermalResistance R_heat = T_ref_100 / Q_el_100
          "Resistance to heat flow";
        constant SI.Temperature T_ref_100 = 100;
      public
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort annotation(Placement(transformation(extent = {{54, 52}, {66, 64}}), iconTransformation(extent = {{90, -100}, {110, -80}})));
        Modelica.Blocks.Interfaces.RealInput u_q annotation(Placement(transformation(extent = {{-198, -48}, {-182, -32}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {-100, -102})));
        Components.Converter converter annotation(Placement(transformation(extent = {{-176, -30}, {-156, -10}})));
        Components.MotorInductance inductance(L_D = L_D) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-90, -40})));
        Components.MotorResistor resistance(R_ref = R_ph_20, alpha = 0.0039, T_ref = 293.15) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-90, 40})));
        Components.EMF airgap(n_polePairs = n_polePairs, k_T_0_100 = k_T_0_100, M_0_100 = M_0_100, I_0_100 = I_0_100, M_0_60 = M_0_60, I_0_60 = I_0_60, M_max = M_max, I_max = I_max) annotation(Placement(transformation(extent = {{-70, -30}, {-50, -10}})));
        Components.Friction friction(C_F = C_F) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-34, 10})));
        Modelica.Mechanics.Rotational.Components.Fixed fixed annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-34, 36})));
        RotationalComponents.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent = {{-4, -30}, {16, -10}})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = C, T(start = T_start, fixed = true)) annotation(Placement(transformation(extent = {{-30, 70}, {-10, 90}})));
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermConductor(G = 1 / R_heat) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {36, 58})));
        HelpBlocks.TemperatureSensor tempSensor annotation(Placement(transformation(extent = {{-6, 60}, {14, 80}})));
        Modelica.Mechanics.Rotational.Components.Inertia motorInertia(J = J) annotation(Placement(transformation(extent = {{26, -30}, {46, -10}})));
        RotationalComponents.Sensors.SpeedSensor_rpm wSensor_rpm annotation(Placement(transformation(extent = {{-4, 0}, {16, 20}})));
        Components.MotorCurrentSensor currentSensor annotation(Placement(transformation(extent = {{-130, 30}, {-110, 50}})));
        Components.MotorVoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-110, -20})));
        Modelica.Blocks.Interfaces.RealOutput i_q
          "current in the branch from p to n as output signal"                                         annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-124, 16}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-100, 90})));
        Components.MotorGround motorGround annotation(Placement(transformation(extent = {{-156, -60}, {-136, -40}})));
      initial equation
        util_n_max = 0;
        util_M_max = 0;
        util_T_max = 0;
        util_U_max = 0;
      equation
        //Optimization problem
        con_n_max = S_n * wSensor_rpm.max - n_max;
        con_M_max = S_M * torqueSensor.max - M_max;
        con_T_max = S_T * tempSensor.actual - T_max;
        con_U_max = S_U * voltageSensor.u_link_max - U_eff_max;
        obj_rms_current = currentSensor.i_rms;
        obj_absMax_current = currentSensor.i_max;
        obj_power_average = converter.power_average;
        connect(airgap.flange, torqueSensor.flange_a) annotation(Line(points = {{-50, -20}, {-4, -20}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.flange_b, motorInertia.flange_a) annotation(Line(points = {{16, -20}, {26, -20}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(motorInertia.flange_b, flange) annotation(Line(points = {{46, -20}, {148, -20}, {148, 0}, {190, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(currentSensor.n, resistance.p) annotation(Line(points = {{-110, 40}, {-100, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(friction.support, fixed.flange) annotation(Line(points = {{-34, 20}, {-34, 36}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(friction.flange, torqueSensor.flange_a) annotation(Line(points = {{-34, 0}, {-34, -20}, {-4, -20}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(voltageSensor.n, currentSensor.n) annotation(Line(points = {{-110, -10}, {-110, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(resistance.heatPort, thermConductor.port_b) annotation(Line(points = {{-90, 44}, {-90, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, thermConductor.port_b) annotation(Line(points = {{-20, 70}, {-20, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, tempSensor.port) annotation(Line(points = {{-20, 70}, {-5.4, 70}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(friction.heatPort, thermConductor.port_b) annotation(Line(points = {{-24, 10}, {-20, 10}, {-20, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(wSensor_rpm.flange, torqueSensor.flange_a) annotation(Line(points = {{-4, 10}, {-12, 10}, {-12, -20}, {-4, -20}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(inductance.p, airgap.n) annotation(Line(points = {{-80, -40}, {-60, -40}, {-60, -30}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(airgap.p, resistance.n) annotation(Line(points = {{-60, -10}, {-60, 40}, {-80, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(airgap.heatport, thermConductor.port_b) annotation(Line(points = {{-54, -18}, {-54, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(converter.u_q, u_q) annotation(Line(points = {{-170, -31}, {-170, -40}, {-190, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(heatPort, thermConductor.port_a) annotation(Line(points = {{60, 58}, {46, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(voltageSensor.p, inductance.n) annotation(Line(points = {{-110, -30}, {-110, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(currentSensor.i_q, i_q) annotation(Line(points = {{-122, 30.1}, {-124, 30.1}, {-124, 16}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(motorGround.p, converter.n) annotation(Line(points = {{-146, -40}, {-146, -24.1}, {-156, -24.1}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(motorGround.p, inductance.n) annotation(Line(points = {{-146, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(converter.p, currentSensor.p) annotation(Line(points = {{-156, -16.1}, {-152, -16.1}, {-152, -16}, {-146, -16}, {-146, 39.9}, {-130, 39.9}}, color = {0, 0, 255}, smooth = Smooth.None));
        when terminal() then
          util_n_max = S_n * wSensor_rpm.max / n_max;
          util_M_max = S_M * torqueSensor.max / M_max;
          util_T_max = 1 + con_T_max;
          //this not physically rigerous. Rigerous values require many cylces. This assumes that if temperature decreases in one cycle by 1 K, we have 0 % utilization. But of course this way also negative values are possible
          util_U_max = S_U * voltageSensor.u_link_max / U_eff_max;
        end when;
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-200,
                  -100},{200,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics={  Text(extent = {{-40, 140}, {240, 100}}, lineColor = {0, 0, 255}, textString = "%name"), Ellipse(extent = {{20, 80}, {180, -80}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{60, 40}, {140, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{80, 28}, {120, 8}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{80, 28}, {100, 8}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{60, 0}, {140, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "M"), Text(extent = {{60, -40}, {140, -80}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "3~"), Rectangle(extent = {{-180, 80}, {-20, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-166, 40}, {-34, -40}}, lineColor = {0, 0, 0}, textString = "Converter"), Line(points = {{-20, 40}, {30, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, 0}, {20, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, -40}, {30, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 0}, {-180, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 40}, {-180, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, -40}, {-180, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-340, 44}, {-224, -44}}, lineColor = {0, 0, 0}, textString = "Power
Supply"), Text(extent = {{-250, -80}, {-50, -118}}, lineColor = {0, 0, 0}, textString = "Uq"), Text(extent = {{-242, 124}, {-42, 86}}, lineColor = {0, 0, 0}, textString = "Iq")}));
      end PMSM;

      model LM
        "Model of a linear motor that can be parametrized with typical supplier data such as from Siemens 1FN3"
        parameter SI.Force F_N "Rated force" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_n "Rated current" annotation(Dialog(group = "Engineering data"));
        parameter SI.Mass m "Mass of primary motor part" annotation(Dialog(group = "Engineering data"));
        parameter SI.Force F_max "Maximum torque" annotation(Dialog(group = "Limiting data"));
        parameter SI.Current I_max "Maximum current" annotation(Dialog(group = "Limiting data"));
        parameter SI.Voltage U_eff_max = 380 "Voltage limit of motor module" annotation(Dialog(group = "Limiting data"));
        parameter SI.ElectricalForceConstant k_F_0_20
          "Rated force constant at 20degC"                                             annotation(Dialog(group = "Physical constants"));
        parameter SI.Resistance R_ph_20 "Winding resistance at 20degC" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inductance L_D "Rotating field inductance" annotation(Dialog(group = "Physical constants"));
        parameter SI.Time T_th(displayUnit = "min") "Thermal time constant" annotation(Dialog(group = "Physical constants"));
        parameter SI.Length tau_p "Pole width" annotation(Dialog(group = "Physical constants"));
        parameter SI.Temperature T_start "Start temperature" annotation(Dialog(group = "Start conditions"));
        parameter Real S_F = 1
          "Safety factor (>=1) for torque (S_n * M <= M_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_U = 1
          "Safety factor (>=1) for voltage (S_n * v <= v_max)"                      annotation(Dialog(tab = "Safety factors"));
        Real con_F_max;
        Real con_U_max;
        Real obj_rms_current;
        Real obj_absMax_current;
      protected
        parameter SI.HeatCapacity C = T_th / R_heat "Heat Capacity";
        parameter SI.HeatFlowRate Q_el_N = 3 * I_n ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at M = M_0_100";
        parameter SI.ThermalResistance R_heat = T_ref_100 / Q_el_N
          "Resistance to heat flow";
        constant SI.Temperature T_ref_100 = 100;
      public
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort annotation(Placement(transformation(extent = {{54, 52}, {66, 64}}), iconTransformation(extent = {{90, -100}, {110, -80}})));
        Modelica.Blocks.Interfaces.RealInput u_q annotation(Placement(transformation(extent = {{-198, -48}, {-182, -32}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {-100, -102})));
        Components.Converter converter annotation(Placement(transformation(extent = {{-176, -30}, {-156, -10}})));
        Components.MotorInductance inductance(L_D = L_D) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-90, -40})));
        Components.MotorResistor resistance(R_ref = R_ph_20, alpha = 0.0039, T_ref = 293.15) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-90, 40})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = C, T(start=
                T_start, fixed=true))                                                                   annotation(Placement(transformation(extent = {{-30, 70}, {-10, 90}})));
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermConductor(G = 1 / R_heat) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {36, 58})));
        HelpBlocks.TemperatureSensor tempSensor annotation(Placement(transformation(extent = {{-6, 60}, {14, 80}})));
        Components.MotorCurrentSensor currentSensor annotation(Placement(transformation(extent = {{-130, 30}, {-110, 50}})));
        Components.MotorVoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-110, -20})));
        Modelica.Blocks.Interfaces.RealOutput i_q
          "current in the branch from p to n as output signal"                                         annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-124, 16}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-100, 90})));
        Components.MotorGround motorGround annotation(Placement(transformation(extent = {{-156, -60}, {-136, -40}})));
        Components.EMF_linear emf(k_F_0_20 = k_F_0_20, F_N = F_N, I_N = I_n, F_max = F_max, I_max = I_max, tau_p = tau_p) annotation(Placement(transformation(extent = {{-70, -30}, {-50, -10}})));
        TranslationalComponents.Sensors.TranslationalSensor transSensor annotation(Placement(transformation(extent = {{-28, -30}, {-8, -10}})));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(Placement(transformation(extent = {{10, -30}, {30, -10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b1
          "Right flange of translational component"                                                              annotation(Placement(transformation(extent = {{170, -10}, {190, 10}}), iconTransformation(extent = {{170, -10}, {190, 10}})));
      equation
        //Optimization problem
        con_F_max = S_F * transSensor.f_max - F_max;
        con_U_max = S_U * voltageSensor.u_link_max - U_eff_max;
        obj_rms_current = currentSensor.i_rms;
        obj_absMax_current = currentSensor.i_max;
        connect(currentSensor.n, resistance.p) annotation(Line(points = {{-110, 40}, {-100, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(voltageSensor.n, currentSensor.n) annotation(Line(points = {{-110, -10}, {-110, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(resistance.heatPort, thermConductor.port_b) annotation(Line(points = {{-90, 44}, {-90, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, thermConductor.port_b) annotation(Line(points = {{-20, 70}, {-20, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, tempSensor.port) annotation(Line(points={{-20,70},
                {-5.4,70}},                                                                                                    color = {191, 0, 0}, smooth = Smooth.None));
        connect(converter.u_q, u_q) annotation(Line(points = {{-170, -31}, {-170, -40}, {-190, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(heatPort, thermConductor.port_a) annotation(Line(points = {{60, 58}, {46, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(voltageSensor.p, inductance.n) annotation(Line(points = {{-110, -30}, {-110, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(currentSensor.i_q, i_q) annotation(Line(points = {{-122, 30.1}, {-124, 30.1}, {-124, 16}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(motorGround.p, converter.n) annotation(Line(points = {{-146, -40}, {-146, -24.1}, {-156, -24.1}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(motorGround.p, inductance.n) annotation(Line(points = {{-146, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(converter.p, currentSensor.p) annotation(Line(points = {{-156, -16.1}, {-152, -16.1}, {-152, -16}, {-146, -16}, {-146, 39.9}, {-130, 39.9}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(inductance.p, emf.n) annotation(Line(points = {{-80, -40}, {-60, -40}, {-60, -30}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(emf.p, resistance.n) annotation(Line(points = {{-60, -10}, {-60, 40}, {-80, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(emf.flange, transSensor.flange_a) annotation(Line(points = {{-54, -20}, {-28.2, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(transSensor.flange_b, mass.flange_a) annotation(Line(points = {{-8, -20}, {10, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(mass.flange_b, flange_b1) annotation(Line(points = {{30, -20}, {106, -20}, {106, 0}, {180, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(emf.heatport, thermConductor.port_b) annotation(Line(points = {{-65, -20}, {-65, 58}, {26, 58}}, color = {191, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics={  Ellipse(extent=  {{20, 80}, {180, -80}}, lineColor=  {0, 0, 0}), Text(extent=  {{-40, 140}, {240, 100}}, lineColor=  {0, 0, 255}, textString=  "%name"), Text(extent=  {{40, 30}, {160, -30}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "M"), Rectangle(extent=  {{-180, 80}, {-20, -80}}, lineColor=  {0, 0, 0}), Text(extent=  {{-166, 40}, {-34, -40}}, lineColor=  {0, 0, 0}, textString=  "Converter"), Line(points=  {{-20, 40}, {30, 40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-20, 0}, {20, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-20, -40}, {30, -40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-218, 0}, {-180, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-218, 40}, {-180, 40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-218, -40}, {-180, -40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Text(extent=  {{-340, 44}, {-224, -44}}, lineColor=  {0, 0, 0}, textString=  "Power
Supply"), Text(extent=  {{-250, -80}, {-50, -118}}, lineColor=  {0, 0, 0}, textString=  "Uq"), Text(extent=  {{-242, 124}, {-42, 86}}, lineColor=  {0, 0, 0}, textString=  "Iq"), Line(points=  {{20, -60}, {182, -60}}, color=  {0, 0, 0}, smooth=  Smooth.None, thickness=  0.5)}));
      end LM;

      model LM_cooling "Model of a linear motor with cooling"
        parameter SI.Force F_N "Rated force" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_n "Rated current" annotation(Dialog(group = "Engineering data"));
        parameter SI.Mass m "Mass of primary motor part" annotation(Dialog(group = "Engineering data"));
        parameter SI.Force F_max "Maximum torque" annotation(Dialog(group = "Limiting data"));
        parameter SI.Current I_max "Maximum current" annotation(Dialog(group = "Limiting data"));
        parameter SI.Voltage U_eff_max = 380 "Voltage limit of motor module" annotation(Dialog(group = "Limiting data"));
        parameter SI.ElectricalForceConstant k_F_0_20
          "Rated force constant at 20degC"                                             annotation(Dialog(group = "Physical constants"));
        parameter SI.Resistance R_ph_20 "Winding resistance at 20degC" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inductance L_D "Rotating field inductance" annotation(Dialog(group = "Physical constants"));
        parameter SI.Time T_th(displayUnit = "min") "Thermal time constant" annotation(Dialog(group = "Physical constants"));
        parameter SI.Length tau_p "Pole width" annotation(Dialog(group = "Physical constants"));
        parameter SI.Temperature T_start = 393.15 "Start temperature"  annotation(Dialog(group = "Start conditions"));
        parameter SI.Temperature T_max(displayUnit = "degC") = 393.15
          "Permissable temperature"                                                             annotation(Dialog(group = "Limiting data"));
        parameter Cooling.Units.FlowResistanceL R_flow = Delta_p_N / Q_N ^ 1.75
          "Resistance to flow of pipe"                                                                       annotation(Dialog(group = "Cooling system"));
        parameter SI.VolumeFlowRate Q_N(displayUnit = "l/min")
          "Nominal volume flow"                                                      annotation(Dialog(group = "Cooling system"));
        parameter SI.Pressure Delta_p_N(displayUnit = "bar")
          "Nominal pressure loss"                                                    annotation(Dialog(group = "Cooling system"));
        replaceable model FluidProp = Cooling.Interfaces.Water
          "constants of fluid models"                                                      annotation(Dialog(group = "Cooling system"));
        FluidProp fluidProp;
        parameter Real S_F = 1
          "Safety factor (>=1) for torque (S_n * M <= M_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_U = 1
          "Safety factor (>=1) for voltage (S_n * v <= v_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_T = 1
          "Safety factor (>=1) for temperature (S_T * T <= T_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter SI.HeatFlowRate Q_el_N = 3 * I_n ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at M = M_0_100";
        Real con_F_max;
        Real con_U_max;
        Real con_T_max;
        Real util_F_max;
        Real util_T_max;
        Real util_U_max;
        Real obj_rms_current;
        Real obj_absMax_current;
        parameter SI.ThermalResistance R_therm = -1 / (Q_N * fluidProp.rho * fluidProp.c_p * log(1 - Q_el_N / (Q_N * fluidProp.rho * fluidProp.c_p * T_ref_100)))
          "Thermal resistance";
        Modelica.Blocks.Interfaces.RealInput u_q annotation(Placement(transformation(extent = {{-198, -48}, {-182, -32}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {-100, -102})));
        Components.Converter converter annotation(Placement(transformation(extent = {{-176, -30}, {-156, -10}})));
        Components.MotorInductance inductance(L_D = L_D) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-90, -40})));
        Components.MotorResistor resistance(R_ref = R_ph_20, alpha = 0.0039, T_ref = 293.15) annotation(Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = 0, origin = {-90, 40})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = C, T(start=
                T_start, fixed=true))                                                                   annotation(Placement(transformation(extent = {{-66, 42}, {-46, 62}})));
        HelpBlocks.TemperatureSensor tempSensor annotation(Placement(transformation(extent = {{-16, 16}, {4, 36}})));
        Components.MotorCurrentSensor currentSensor annotation(Placement(transformation(extent = {{-130, 30}, {-110, 50}})));
        Components.MotorVoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-110, -20})));
        Modelica.Blocks.Interfaces.RealOutput i_q
          "current in the branch from p to n as output signal"                                         annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-124, 16}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-100, 90})));
        Components.MotorGround motorGround annotation(Placement(transformation(extent = {{-156, -60}, {-136, -40}})));
        Components.EMF_linear emf_linear(k_F_0_20 = k_F_0_20, F_N = F_N, I_N = I_n, F_max = F_max, I_max = I_max, tau_p = tau_p) annotation(Placement(transformation(extent = {{-70, -30}, {-50, -10}})));
        TranslationalComponents.Sensors.TranslationalSensor transSensor annotation(Placement(transformation(extent = {{-40, -30}, {-20, -10}})));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(Placement(transformation(extent = {{-10, -30}, {10, -10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b f
          "Right flange of translational component"                                                      annotation(Placement(transformation(extent = {{12, -26}, {24, -14}}), iconTransformation(extent = {{180, -10}, {200, 10}})));
        Cooling.CoolingChannels coolingChannels(R_flow = R_flow, R_therm = R_therm) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-22, 50})));
        Cooling.Interfaces.HydraulicPort_b coolingPort_b annotation(Placement(transformation(extent = {{154, 54}, {168, 68}}), iconTransformation(extent = {{148, 48}, {168, 68}})));
        Cooling.Interfaces.HydraulicPort_a coolingPort_a annotation(Placement(transformation(extent = {{38, 54}, {52, 68}}), iconTransformation(extent = {{32, 48}, {52, 68}})));
        Modelica.Blocks.Interfaces.RealOutput actual1 annotation(Placement(transformation(extent = {{96, 80}, {116, 100}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {100, 88})));
      protected
        parameter SI.HeatCapacity C = T_th / R_therm "Heat Capacity";
        constant SI.Temperature T_ref_100 = 100;
      initial equation
        util_F_max = 0;
        util_T_max = 0;
        util_U_max = 0;
      equation
        //Optimization problem
        con_F_max = S_F * transSensor.f_max - F_max;
        con_U_max = S_U * voltageSensor.u_link_max - U_eff_max;
        con_T_max = S_T * tempSensor.actual - T_max;
        obj_rms_current = currentSensor.i_rms;
        obj_absMax_current = currentSensor.i_max;
        when terminal() then
          util_F_max = S_F * transSensor.f_max / F_max;
          util_T_max = 1 + con_T_max;
          //this not physically rigerous. Rigerous values require many cylces. This assumes that if temperature decreases in one cycle by 1 K, we have 0 % utilization. But of course this way also negative values are possible
          util_U_max = S_U * voltageSensor.u_link_max / U_eff_max;
        end when;
        connect(currentSensor.n, resistance.p) annotation(Line(points = {{-110, 40}, {-100, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(voltageSensor.n, currentSensor.n) annotation(Line(points = {{-110, -10}, {-110, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(converter.u_q, u_q) annotation(Line(points = {{-170, -31}, {-170, -40}, {-190, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(voltageSensor.p, inductance.n) annotation(Line(points = {{-110, -30}, {-110, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(currentSensor.i_q, i_q) annotation(Line(points = {{-122, 30.1}, {-124, 30.1}, {-124, 16}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(motorGround.p, converter.n) annotation(Line(points = {{-146, -40}, {-146, -24.1}, {-156, -24.1}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(motorGround.p, inductance.n) annotation(Line(points = {{-146, -40}, {-100, -40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(converter.p, currentSensor.p) annotation(Line(points = {{-156, -16.1}, {-152, -16.1}, {-152, -16}, {-146, -16}, {-146, 39.9}, {-130, 39.9}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(inductance.p, emf_linear.n) annotation(Line(points = {{-80, -40}, {-60, -40}, {-60, -30}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(emf_linear.p, resistance.n) annotation(Line(points = {{-60, -10}, {-60, 40}, {-80, 40}}, color = {0, 0, 255}, smooth = Smooth.None));
        connect(emf_linear.flange, transSensor.flange_a) annotation(Line(points = {{-54, -20}, {-40.2, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(emf_linear.heatport, resistance.heatPort) annotation(Line(points = {{-65, -20}, {-90, -20}, {-90, 36}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, resistance.heatPort) annotation(Line(points = {{-56, 42}, {-56, 26}, {-90, 26}, {-90, 36}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(tempSensor.port, resistance.heatPort) annotation(Line(points = {{-15.4, 26}, {-90, 26}, {-90, 36}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(coolingChannels.heatport_a, resistance.heatPort) annotation(Line(points = {{-22, 39}, {-22, 26}, {-90, 26}, {-90, 36}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(coolingChannels.port_a, coolingPort_a) annotation(Line(points = {{-32.4, 50}, {45, 50}, {45, 61}}, color = {0, 0, 0}, pattern = LinePattern.Solid, thickness = 0.25, smooth = Smooth.None));
        connect(coolingChannels.port_b, coolingPort_b) annotation(Line(points = {{-11.6, 50}, {161, 50}, {161, 61}}, color = {0, 0, 0}, pattern = LinePattern.Solid, smooth = Smooth.None));
        connect(mass.flange_a, transSensor.flange_b) annotation(Line(points = {{-10, -20}, {-20, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(mass.flange_b, f) annotation(Line(points = {{10, -20}, {18, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(tempSensor.actual, actual1) annotation(Line(points = {{4.6, 27}, {26.3, 27}, {26.3, 90}, {106, 90}}, color = {0, 0, 127}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}), graphics={  Ellipse(extent = {{20, 80}, {180, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-120, -100}, {340, -140}}, lineColor = {0, 0, 255}, textString = "%name"), Text(extent = {{40, 30}, {160, -30}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "M"), Rectangle(extent = {{-180, 80}, {-20, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-166, 40}, {-34, -40}}, lineColor = {0, 0, 0}, textString = "Converter"), Line(points = {{-20, 40}, {30, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, 0}, {20, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, -40}, {30, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 0}, {-180, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 40}, {-180, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, -40}, {-180, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-340, 44}, {-224, -44}}, lineColor = {0, 0, 0}, textString = "Power
Supply"), Text(extent = {{-250, -80}, {-50, -118}}, lineColor = {0, 0, 0}, textString = "Uq"), Text(extent = {{-242, 124}, {-42, 86}}, lineColor = {0, 0, 0}, textString = "Iq"), Line(points = {{20, -60}, {182, -60}}, color = {0, 0, 0}, smooth = Smooth.None, thickness = 0.5), Text(extent = {{40, 82}, {160, 40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "T")}));
      end LM_cooling;
    end Motors;

    package Motors_withAutoTunedCurrentController
      model PMSM "Simple Model of a synchronous servo motor"
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_n
          "Rated speed"                                                           annotation(Dialog(group = "Engineering data"));
        parameter Integer n_polePairs "Number of pole pairs" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_n "Rated torque (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_n "Rated current (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_0_100 "Stall torque (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_0_100 "Stall current (100 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_0_60 "Stall torque (60 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Current I_0_60 "Stall current (60 K)" annotation(Dialog(group = "Engineering data"));
        parameter SI.Inertia J "Motor inertia" annotation(Dialog(group = "Engineering data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_max
          "Max. permissible speed (mech.)"                                                             annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_max "Maximum torque" annotation(Dialog(group = "Limiting data"));
        parameter SI.Current I_max "Maximum current" annotation(Dialog(group = "Limiting data"));
        parameter SI.Temperature T_max(displayUnit = "degC") = 413.15
          "Permissable temperature"                                                             annotation(Dialog(group = "Limiting data"));
        parameter SI.Voltage U_eff_max = 380 "Voltage limit of motor module" annotation(Dialog(group = "Limiting data"));
        parameter SI.ElectricalTorqueConstant k_T_0_100 = M_0_100 / I_0_100
          "Rated torque constant"                                                                   annotation(Dialog(group = "Physical constants"));
        parameter SI.Resistance R_ph_20 "Winding resistance at 20degC" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inductance L_D "Rotating field inductance" annotation(Dialog(group = "Physical constants"));
        parameter SI.Time T_th(displayUnit = "min") "Thermal time constant" annotation(Dialog(group = "Physical constants"));
        parameter Real C_F = (T_ref_100 / R_heat - Q_R_N) / SI.Conversions.from_rpm(n_n) ^ 1.5
          "Constant relating friction torque and square root of angular velocity"
                                                                                                              annotation(Dialog(group = "Physical constants"));
        parameter SI.Time T_AT(displayUnit = "ms") = 125e-6
          "Sampling time for current control loop"                                                    annotation(Dialog(group = "Physical constants"));
        parameter SI.Temperature T_start "Start temperature" annotation(Dialog(group = "Start conditions"));
        parameter Real S_n = 1
          "Safety factor (>=1) for rotation speed (S_n * n <= n_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_M = 1
          "Safety factor (>=1) for torque (S_n * M <= M_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_U = 1
          "Safety factor (>=1) for voltage (S_n * u <= u_max)"                      annotation(Dialog(tab = "Safety factors"));
        parameter Real S_T = 1
          "Safety factor (>=1) for temperature (S_n * T <= T_max)"                      annotation(Dialog(tab = "Safety factors"));
        Real con_n_max;
        Real con_M_max;
        Real con_T_max "End of period value, for simulation of one period";
        Real con_U_max;
        Real util_n_max;
        Real util_M_max;
        Real util_T_max;
        Real util_U_max;
        Real obj_rms_current;
        Real obj_absMax_current;
        Real obj_power_average;
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange
          "Flange of shaft"                                                        annotation(Placement(transformation(extent = {{180, -10}, {200, 10}}), iconTransformation(extent = {{180, -10}, {200, 10}})));
      protected
        parameter SI.HeatCapacity C = T_th / R_heat "Heat Capacity";
        parameter SI.HeatFlowRate Q_el_100 = 3 * I_0_100 ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at M = M_0_100";
        parameter SI.HeatFlowRate Q_R_N = 3 * (M_n / k_T_0_100) ^ 2 * (1 + 100 * 0.0039) * R_ph_20
          "Electrical heat at nominal point";
        parameter SI.ThermalResistance R_heat = T_ref_100 / Q_el_100
          "Resistance to heat flow";
        constant SI.Temperature T_ref_100 = 100;
      public
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort annotation(Placement(transformation(extent={{130,62},
                  {142,74}}),                                                                                                    iconTransformation(extent = {{90, -100}, {110, -80}})));
        Components.Converter converter annotation(Placement(transformation(extent={{-100,
                  -20},{-80,0}})));
        Components.MotorInductance inductance(L_D = L_D) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={-14,-30})));
        Components.MotorResistor resistance(R_ref = R_ph_20, alpha = 0.0039, T_ref = 293.15) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin={-14,50})));
        Components.EMF airgap(n_polePairs = n_polePairs, k_T_0_100 = k_T_0_100, M_0_100 = M_0_100, I_0_100 = I_0_100, M_0_60 = M_0_60, I_0_60 = I_0_60, M_max = M_max, I_max = I_max) annotation(Placement(transformation(extent={{6,-20},
                  {26,0}})));
        Components.Friction friction(C_F = C_F) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={42,20})));
        Modelica.Mechanics.Rotational.Components.Fixed fixed annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={42,46})));
        RotationalComponents.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent={{72,-20},
                  {92,0}})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = C, T(start = T_start, fixed = true)) annotation(Placement(transformation(extent={{46,80},
                  {66,100}})));
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermConductor(G = 1 / R_heat) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={112,68})));
        HelpBlocks.TemperatureSensor tempSensor annotation(Placement(transformation(extent={{70,70},
                  {90,90}})));
        Modelica.Mechanics.Rotational.Components.Inertia motorInertia(J = J) annotation(Placement(transformation(extent={{102,-20},
                  {122,0}})));
        RotationalComponents.Sensors.SpeedSensor_rpm wSensor_rpm annotation(Placement(transformation(extent={{72,10},
                  {92,30}})));
        Components.MotorCurrentSensor currentSensor annotation(Placement(transformation(extent={{-54,40},
                  {-34,60}})));
        Components.MotorVoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin={-34,-10})));
        Components.MotorGround motorGround annotation(Placement(transformation(extent={{-76,-50},
                  {-56,-30}})));
        Modelica.Blocks.Continuous.PI PI_I(                       T=L_D/R_ph_20, k=L_D/(2
              *1.5*T_AT))
          "Control according to amplitude optimum. k = Ts/(2*Ks*T_omega), where s denotes the controlled system with Ks = 1/R and T_omega = 1.5 * T_AT (sampling time)"
                                                                  annotation(Placement(transformation(extent={{-120,
                  -80},{-100,-60}})));
        HelpBlocks.Feedback_mirror I_feedback annotation(Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 180, origin={-142,-38})));
        Modelica.Blocks.Interfaces.RealInput i_set
          annotation (Placement(transformation(extent={{-124,-114},{-84,-74}}),
              iconTransformation(
              extent={{-20,-20},{20,20}},
              rotation=90,
              origin={-100,-100})));
        Components.voltageLimit voltageLimit(U_eff_max=U_eff_max)
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-94,-44})));
        Modelica.Blocks.Nonlinear.Limiter limiter(uMax=I_max)
          annotation (Placement(transformation(extent={{-184,10},{-164,30}})));
        Modelica.Blocks.Continuous.FirstOrder T_GI(k=1, T=0.08*L_D/R_ph_20)
          "0,08 is just a guess value to avoid limits"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-160,-14})));
        Modelica.Blocks.Math.Gain gain(k=I_0_100/M_0_100) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-194,-60})));
      initial equation
        util_n_max = 0;
        util_M_max = 0;
        util_T_max = 0;
        util_U_max = 0;
      equation
        //Optimization problem
        con_n_max = S_n * wSensor_rpm.max - n_max;
        con_M_max = S_M * torqueSensor.max - M_max;
        con_T_max = S_T * tempSensor.actual - T_max;
        con_U_max = S_U * voltageSensor.u_link_max - U_eff_max;
        obj_rms_current = currentSensor.i_rms;
        obj_absMax_current = currentSensor.i_max;
        obj_power_average = converter.power_average;
        connect(airgap.flange, torqueSensor.flange_a) annotation(Line(points={{26,-10},
                {72,-10}},                                                                              color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.flange_b, motorInertia.flange_a) annotation(Line(points={{92,-10},
                {102,-10}},                                                                                    color = {0, 0, 0}, smooth = Smooth.None));
        connect(motorInertia.flange_b, flange) annotation(Line(points={{122,-10},{178,
                -10},{178,0},{190,0}},                                                                               color = {0, 0, 0}, smooth = Smooth.None));
        connect(currentSensor.n, resistance.p) annotation(Line(points={{-34,50},{-24,50}},        color = {0, 0, 255}, smooth = Smooth.None));
        connect(friction.support, fixed.flange) annotation(Line(points={{42,30},{42,46}},        color = {0, 0, 0}, smooth = Smooth.None));
        connect(friction.flange, torqueSensor.flange_a) annotation(Line(points={{42,10},
                {42,-10},{72,-10}},                                                                                 color = {0, 0, 0}, smooth = Smooth.None));
        connect(voltageSensor.n, currentSensor.n) annotation(Line(points={{-34,0},{-34,
                50}},                                                                                 color = {0, 0, 255}, smooth = Smooth.None));
        connect(resistance.heatPort, thermConductor.port_b) annotation(Line(points={{-14,54},
                {-14,68},{102,68}},                                                                                    color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, thermConductor.port_b) annotation(Line(points={{56,80},
                {56,68},{102,68}},                                                                                    color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor.port, tempSensor.port) annotation(Line(points={{56,80},{
                70.6,80}},                                                                             color = {191, 0, 0}, smooth = Smooth.None));
        connect(friction.heatPort, thermConductor.port_b) annotation(Line(points={{52,20},
                {56,20},{56,68},{102,68}},                                                                                      color = {191, 0, 0}, smooth = Smooth.None));
        connect(wSensor_rpm.flange, torqueSensor.flange_a) annotation(Line(points={{72,20},
                {64,20},{64,-10},{72,-10}},                                                                                       color = {0, 0, 0}, smooth = Smooth.None));
        connect(inductance.p, airgap.n) annotation(Line(points={{-4,-30},{16,-30},{16,
                -20}},                                                                                 color = {0, 0, 255}, smooth = Smooth.None));
        connect(airgap.p, resistance.n) annotation(Line(points={{16,0},{16,50},{-4,50}},             color = {0, 0, 255}, smooth = Smooth.None));
        connect(airgap.heatport, thermConductor.port_b) annotation(Line(points={{22,-8},
                {22,68},{102,68}},                                                                                  color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatPort, thermConductor.port_a) annotation(Line(points={{136,68},{122,
                68}},                                                                           color = {191, 0, 0}, smooth = Smooth.None));
        connect(voltageSensor.p, inductance.n) annotation(Line(points={{-34,-20},{-34,
                -30},{-24,-30}},                                                                                 color = {0, 0, 255}, smooth = Smooth.None));
        connect(motorGround.p, converter.n) annotation(Line(points={{-66,-30},{
                -66,-14.1},{-80,-14.1}},                                                                          color = {0, 0, 255}, smooth = Smooth.None));
        connect(motorGround.p, inductance.n) annotation(Line(points={{-66,-30},
                {-24,-30}},                                                                       color = {0, 0, 255}, smooth = Smooth.None));
        connect(converter.p, currentSensor.p) annotation(Line(points={{-80,-6.1},
                {-76,-6.1},{-76,-6},{-70,-6},{-70,49.9},{-54,49.9}},                                                                                         color = {0, 0, 255}, smooth = Smooth.None));
        when terminal() then
          util_n_max = S_n * wSensor_rpm.max / n_max;
          util_M_max = S_M * torqueSensor.max / M_max;
          util_T_max = 1 + con_T_max;
          //this not physically rigerous. Rigerous values require many cylces. This assumes that if temperature decreases in one cycle by 1 K, we have 0 % utilization. But of course this way also negative values are possible
          util_U_max = S_U * voltageSensor.u_link_max / U_eff_max;
        end when;
        connect(PI_I.u, I_feedback.y) annotation (Line(
            points={{-122,-70},{-128,-70},{-128,-38},{-133,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(I_feedback.u2, currentSensor.i_q) annotation (Line(
            points={{-142,-30},{-142,30},{-46,30},{-46,40.1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageLimit.y, converter.u_q) annotation (Line(
            points={{-94,-33},{-94,-21}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageLimit.u_d, voltageSensor.u_d_out) annotation (Line(
            points={{-86,-56},{-86,-74},{-16,-74},{-16,-13},{-24,-13}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(PI_I.y, voltageLimit.u) annotation (Line(
            points={{-99,-70},{-94,-70},{-94,-56}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T_GI.y, I_feedback.u1) annotation (Line(
            points={{-160,-25},{-160,-38},{-150,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T_GI.u, limiter.y) annotation (Line(
            points={{-160,-2},{-160,20},{-163,20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(gain.y, limiter.u) annotation (Line(
            points={{-194,-49},{-194,20},{-186,20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(gain.u, i_set) annotation (Line(
            points={{-194,-72},{-194,-94},{-104,-94}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-200,
                  -100},{200,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio=false,   extent={{-200,
                  -100},{200,100}}),                                                                                                    graphics={  Text(extent = {{-40, 140}, {240, 100}}, lineColor = {0, 0, 255}, textString = "%name"), Ellipse(extent = {{20, 80}, {180, -80}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{60, 40}, {140, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{80, 28}, {120, 8}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{80, 28}, {100, 8}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{60, 0}, {140, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "M"), Text(extent = {{60, -40}, {140, -80}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid, textString = "3~"), Rectangle(extent = {{-180, 80}, {-20, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-166, 40}, {-34, -40}}, lineColor = {0, 0, 0}, textString = "Converter"), Line(points = {{-20, 40}, {30, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, 0}, {20, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-20, -40}, {30, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 0}, {-180, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, 40}, {-180, 40}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-218, -40}, {-180, -40}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-340, 44}, {-224, -44}}, lineColor = {0, 0, 0}, textString = "Power
Supply"), Text(extent = {{-250, -80}, {-50, -118}}, lineColor={0,0,0},
                textString="Iq")}));
      end PMSM;

    end Motors_withAutoTunedCurrentController;

    package Components
      model HeatTransfer "Comination of capacitor and thermal resistance"
        parameter SI.HeatCapacity C "Heat Capacity";
        parameter SI.ThermalResistance R_heat "Resistance to heat flow";
        parameter SI.Temperature T_start "Start temperature";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(transformation(extent = {{-120, -10}, {-100, 10}}), iconTransformation(extent = {{-120, -10}, {-100, 10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1(T(start = T_start), C = C) annotation(Placement(transformation(extent = {{-12, 18}, {8, 38}})));
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G = 1 / R_heat) annotation(Placement(transformation(extent = {{22, -10}, {42, 10}})));
      equation
        connect(heatCapacitor1.port, thermalConductor1.port_a) annotation(Line(points = {{-2, 18}, {-2, 0}, {22, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(thermalConductor1.port_b, port_b) annotation(Line(points = {{42, 0}, {110, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor1.port, port_a) annotation(Line(points = {{-2, 18}, {-2, 0}, {-110, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-88, 70}, {92, -70}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.None, fillColor=  {192, 192, 192},
                  fillPattern=                                                                                                    FillPattern.Backward), Line(points=  {{-88, 70}, {-88, -70}}, color=  {0, 0, 0}, thickness=  0.5), Line(points=  {{92, 70}, {92, -70}}, color=  {0, 0, 0}, thickness=  0.5), Text(extent=  {{-148, 115}, {152, 75}}, textString=  "%name", lineColor=  {0, 0, 255})}));
      end HeatTransfer;

      model HeatTransfer_2C "Model with two heat capcitors"
        parameter SI.HeatCapacity C "Heat Capacity";
        parameter SI.ThermalResistance R_heat "Resistance to heat flow";
        parameter SI.Temperature T_start "Start temperature";
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 2 * 1 / R_heat) annotation(Placement(transformation(extent = {{-40, -10}, {-20, 10}})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(T(start = T_start), C = 0.2 * C) annotation(Placement(transformation(extent = {{-60, 16}, {-40, 36}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(transformation(extent = {{-120, -10}, {-100, 10}}), iconTransformation(extent = {{-120, -10}, {-100, 10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1(T(start = T_start), C = 0.8 * C) annotation(Placement(transformation(extent = {{-12, 18}, {8, 38}})));
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G = 2 * 1 / R_heat) annotation(Placement(transformation(extent = {{22, -10}, {42, 10}})));
      equation
        connect(heatCapacitor.port, thermalConductor.port_a) annotation(Line(points = {{-50, 16}, {-50, 0}, {-40, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(thermalConductor.port_a, port_a) annotation(Line(points = {{-40, 0}, {-110, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(heatCapacitor1.port, thermalConductor1.port_a) annotation(Line(points = {{-2, 18}, {-2, 0}, {22, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(thermalConductor.port_b, thermalConductor1.port_a) annotation(Line(points = {{-20, 0}, {22, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        connect(thermalConductor1.port_b, port_b) annotation(Line(points = {{42, 0}, {110, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-88, 70}, {92, -70}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.None, fillColor=  {192, 192, 192},
                  fillPattern=                                                                                                    FillPattern.Backward), Line(points=  {{-88, 70}, {-88, -70}}, color=  {0, 0, 0}, thickness=  0.5), Line(points=  {{92, 70}, {92, -70}}, color=  {0, 0, 0}, thickness=  0.5), Text(extent=  {{-148, 115}, {152, 75}}, textString=  "%name", lineColor=  {0, 0, 255})}));
      end HeatTransfer_2C;

      model Converter
        "Model for a voltage-driven converter. The field producing current i_d is set to zero"

        Modelica.Blocks.Interfaces.RealInput u_q annotation(Placement(transformation(extent = {{-140, 60}, {-100, 100}}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {-80, -220})));
        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{190, 68}, {210, 88}}), iconTransformation(extent = {{190, 68}, {210, 88}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{190, -92}, {210, -72}}), iconTransformation(extent = {{190, -92}, {210, -72}})));
        SI.Voltage u_d;
        SI.Power power "Acutual power";
        SI.Power power_average "Average power";
        SI.Power power_no_recovery_average "Average power if without recovery";
        SI.Energy energy;
        SI.Energy energy_no_recovery;
      protected
        parameter Real startTime(fixed = false);
      initial equation
        startTime = time;
        energy = 0;
        energy_no_recovery = 0;
      equation
        u_d = p.u_d - n.u_d;
        u_q = p.u_q - n.u_q;
        0 = p.i_d + n.i_d;
        0 = p.i_q + n.i_q;
        //0 = p.w - n.w;
        p.i_d = 0; //This is a simplification, as i.d is usually set to zero by current control
        power = -3 * (p.i_q * u_q + p.i_d * u_d);
        der(energy) = power;
        der(energy_no_recovery) = max(power, 0);
        if time <= startTime then
          power_average = 0;
          power_no_recovery_average = 0;
        else
          power_average = energy / (time - startTime);
          power_no_recovery_average = energy_no_recovery / (time - startTime);
        end if;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -200}, {200, 200}}), graphics={  Rectangle(extent=  {{-200, 200}, {200, -200}}, lineColor=  {0, 0, 0}),
                                                                                                    Line(points=  {{200, 200}, {-200, -200}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-160, 100}, {-148.7, 134.2}, {-141.5, 153.1}, {-135.1, 166.4}, {-129.4, 174.6}, {-123.8, 179.1}, {-118.2, 179.8}, {-112.6, 176.6}, {-106.9, 169.7}, {-101.3, 159.4}, {-94.9, 144.1}, {-86.83, 121.2}, {-69.9, 69.2}, {-62.7, 49.8}, {-56.3, 35.8}, {-50.7, 26.9}, {-45, 21.6}, {-39.4, 20}, {-33.8, 22.4}, {-28.1, 28.5}, {-22.5, 38.1}, {-16.1, 52.8}, {-8, 75.2}, {0, 100}}, color=  {0, 0, 0}), Line(points=  {{0, -100}, {11.3, -65.8}, {18.5, -46.9}, {24.9, -33.6}, {30.6, -25.4}, {36.2, -20.9}, {41.8, -20.2}, {47.4, -23.4}, {53.1, -30.3}, {58.7, -40.6}, {65.1, -55.9}, {73.17, -78.8}, {90.1, -130.8}, {97.3, -150.2}, {103.7, -164.2}, {109.3, -173.1}, {115, -178.4}, {120.6, -180}, {126.2, -177.6}, {131.9, -171.5}, {137.5, -161.9}, {143.9, -147.2}, {152, -124.8}, {160, -100}}, color=  {0, 0, 0}), Text(extent=  {{-304, -220}, {-120, -300}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "Uq"), Text(extent=  {{-258, 302}, {260, 200}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{-238, 40}, {-200, 40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-238, 0}, {-200, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-238, -40}, {-200, -40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Text(extent=  {{-698, 64}, {-180, -38}}, lineColor=  {0, 0, 0}, textString=  "power 
supply")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -200}, {200, 200}}), graphics));
      end Converter;

      model EMF "Electromotoric force (electric/mechanic transformer)"

        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Interfaces.Support support if useSupport
          "Support/housing of emf shaft"                                                                      annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{-10, -110}, {10, -90}}), iconTransformation(extent = {{-10, -110}, {10, -90}})));
        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{-10, 90}, {10, 110}}), iconTransformation(extent = {{-10, 90}, {10, 110}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport annotation(Placement(transformation(extent = {{50, 10}, {70, 30}}), iconTransformation(extent = {{50, 10}, {70, 30}})));
        parameter Boolean useSupport = false
          "= true, if support flange enabled, otherwise implicitly grounded"                                    annotation(Evaluate = true, HideResult = true, choices(__Dymola_checkBox = true));
        parameter Integer n_polePairs "Number of pole pairs";
        parameter SI.ElectricalTorqueConstant k_T_0_100
          "Rated torque constant at 100 K overtemperature";
        parameter SI.Torque M_0_100 "Stall torque (100 K)";
        parameter SI.Current I_0_100 "Stall current (100 K)";
        parameter SI.Torque M_0_60 "Stall torque (100 K)";
        parameter SI.Current I_0_60 "Stall current (100 K)";
        parameter SI.Torque M_max "Maximum permissible torque";
        parameter SI.Current I_max "Maximum permissible current";
        SI.Angle phi
          "Angle of shaft flange with respect to support (= flange.phi - support.phi)";
        SI.AngularVelocity w_mech
          "Angular velocity of flange relative to support";
        SI.Current i_q(start = 0) "Current";
        SI.ElectricalTorqueConstant k_T(start = k_T_0_100)
          "Actual torque constant";
        SI.Torque tau "Torque at flange";
        SI.Torque tau_abs "Absolute torque at flange";
        SI.Voltage u_emf "Voltage drop at emf";
        SI.Temperature T "Temperature of motor";
      protected
        parameter Real Delta_k_T_Temp = (M_0_60 / I_0_60 - M_0_100 / I_0_100) / 40
          "Decrease of torque constant with temperature";
        parameter Real Delta_k_T_Torque = (M_0_100 / I_0_100 - M_max / I_max) / (M_max - M_0_100)
          "Decrease of torque constant with torque";
        Real c_Torque(start = 1)
          "Decrease factor of torque constant with torque";
        Real c_Temp(start = 1)
          "Decrease factor of torque constant with temperature";
        Modelica.Mechanics.Rotational.Components.Fixed fixed if not useSupport annotation(Placement(transformation(extent = {{-90, -20}, {-70, 0}})));
        Modelica.Mechanics.Rotational.Interfaces.InternalSupport internalSupport(tau = -flange.tau) annotation(Placement(transformation(extent = {{-90, -10}, {-70, 10}})));
      equation
        heatport.Q_flow = 0; // no heat generation
        heatport.T = T;
        u_emf = p.u_q - n.u_q;
        0 = p.u_d - n.u_d;
        0 = p.i_q + n.i_q;
        0 = p.i_d + n.i_d;
        i_q = p.i_q;
        w_mech * n_polePairs = p.w;
        0 = p.w - n.w;
        phi = flange.phi - internalSupport.phi;
        w_mech = der(phi);
        flange.tau = tau;
        tau_abs = abs(tau);
        k_T * w_mech = 3 * u_emf;
        tau = -k_T * i_q;

        //The torque constant k_T is dependent on temperature and torque
        k_T_0_100 * c_Torque = if tau_abs > M_max then M_max / I_max elseif tau_abs < M_0_100 then k_T_0_100 else k_T_0_100 - Delta_k_T_Torque * (tau_abs - M_0_100);
        //100 K and 0 K overtemperature
        k_T_0_100 * c_Temp = if T > 393.15 then k_T_0_100 elseif T < 293.15 then k_T_0_100 + Delta_k_T_Temp * 100 else k_T_0_100 + Delta_k_T_Temp * (393.15 - T);
        k_T = k_T_0_100 * c_Temp * c_Torque;
        connect(internalSupport.flange, support) annotation(Line(points = {{-80, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(internalSupport.flange, fixed.flange) annotation(Line(points = {{-80, 0}, {-80, -10}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(defaultComponentName = "emf", Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Rectangle(extent=  {{-85, 10}, {-36, -10}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {192, 192, 192}), Line(points=  {{0, 90}, {0, 40}}, color=  {0, 0, 255}), Rectangle(extent=  {{35, 10}, {100, -10}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {192, 192, 192}), Ellipse(extent=  {{-40, 40}, {40, -40}}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, lineColor=  {0, 0, 255}), Line(points=  {{0, -90}, {0, -40}}, color=  {0, 0, 255}), Text(extent=  {{0, -50}, {199, -90}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(visible=  not useSupport, points=  {{-100, -30}, {-40, -30}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-100, -50}, {-80, -30}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-80, -50}, {-60, -30}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-60, -50}, {-40, -30}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-70, -30}, {-70, -10}}, color=  {0, 0, 0})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points=  {{-17, 95}, {-20, 85}, {-23, 95}, {-17, 95}}, lineColor=  {160, 160, 164}, fillColor=  {160, 160, 164},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-20, 110}, {-20, 85}}, color=  {160, 160, 164}), Text(extent=  {{-40, 110}, {-30, 90}}, lineColor=  {160, 160, 164}, textString=  "i"), Line(points=  {{9, 75}, {19, 75}}, color=  {192, 192, 192}), Line(points=  {{-20, -110}, {-20, -85}}, color=  {160, 160, 164}), Polygon(points=  {{-17, -100}, {-20, -110}, {-23, -100}, {-17, -100}}, lineColor=  {160, 160, 164}, fillColor=  {160, 160, 164},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-40, -110}, {-30, -90}}, lineColor=  {160, 160, 164}, textString=  "i"), Line(points=  {{8, -79}, {18, -79}}, color=  {192, 192, 192}), Line(points=  {{14, 80}, {14, 70}}, color=  {192, 192, 192})}), Documentation(info = "<html>
<p>EMF transforms electrical energy into rotational mechanical energy. It is used as basic building block of an electrical motor. The mechanical connector flange can be connected to elements of the Modelica.Mechanics.Rotational library. flange.tau is the cut-torque, flange.phi is the angle at the rotational connection.</p>
</html>", revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Martin Otter<br> initially implemented<br>
       </li>
</ul>
</html>"));
      end EMF;

      model EMF_linear "Electromotoric force (electric/mechanic transformer)"

        Modelica.Mechanics.Translational.Interfaces.Flange_b flange annotation(Placement(transformation(extent = {{50, -10}, {70, 10}}, rotation = 0), iconTransformation(extent = {{50, -10}, {70, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Support support if useSupport
          "Support/housing of emf shaft"                                                                         annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{-10, -110}, {10, -90}}), iconTransformation(extent = {{-10, -110}, {10, -90}})));
        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{-10, 90}, {10, 110}}), iconTransformation(extent = {{-10, 90}, {10, 110}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport annotation(Placement(transformation(extent = {{-60, -10}, {-40, 10}}), iconTransformation(extent = {{-60, -10}, {-40, 10}})));
        parameter Boolean useSupport = false
          "= true, if support flange enabled, otherwise implicitly grounded"                                    annotation(Evaluate = true, HideResult = true, choices(__Dymola_checkBox = true));
        parameter SI.ElectricalForceConstant k_F_0_20
          "Rated force constant at 20 ?C";
        parameter SI.Force F_N "Nomina (120 ?C)";
        parameter SI.Current I_N "Nominal current (120 ?C)";
        parameter SI.Force F_max "Maximum permissible force (120 ?C)";
        parameter SI.Current I_max "Maximum permissible current (120 ?C)";
        parameter SI.Length tau_p "Pole width";
        SI.Velocity v_M "Motor velocity";
        SI.Position s
          "Angle of shaft flange with respect to support (= flange.phi - support.phi)";
        SI.Frequency f_M "Motor frequency of the induced voltage";
        SI.Current i_q(start = 0) "Current";
        SI.ElectricalForceConstant k_F(start = k_F_0_20)
          "Actual force constant";
        SI.Force force "Force at flange";
        SI.Force force_abs "Absolute force at flange";
        SI.Voltage u_emf "Voltage drop at emf";
        SI.Temperature T "Temperature of motor";
      protected
        parameter Real Delta_k_F_Force = (F_N / I_N - F_max / I_max) / (F_max - F_N)
          "Decrease of force constant with torque";
        parameter Real Delta_k_F_Temp = (k_F_0_20 - F_N / I_N) / 100
          "Decrease of force constant with froce";
        Real c_Force(start = 1) "Decrease factor of force constant with force";
        Real c_Temp(start = 1)
          "Decrease factor of force constant with temperature";
        Modelica.Mechanics.Translational.Components.Fixed fixed if not useSupport annotation(Placement(transformation(extent = {{-90, -20}, {-70, 0}})));
        Modelica.Mechanics.Translational.Interfaces.InternalSupport internalSupport(f = -flange.f) annotation(Placement(transformation(extent = {{-90, -10}, {-70, 10}})));
      equation
        heatport.Q_flow = 0;
        heatport.T = T;
        u_emf = p.u_q - n.u_q;
        0 = p.u_d - n.u_d;
        0 = p.i_q + n.i_q;
        0 = p.i_d + n.i_d;
        i_q = p.i_q;
        f_M * 2 * Constants.pi = p.w;
        v_M = 2 * tau_p * f_M;
        //factor 2 (see GROS06, 257ff)
        0 = p.w - n.w;
        s = flange.s - internalSupport.s;
        v_M = der(s);
        flange.f = force;
        force_abs = abs(force);
        k_F * v_M = 3 * u_emf;
        force = -k_F * i_q;
        //The force constant k_T is dependent on temperature and force
        F_N / I_N * c_Force = if force_abs > F_max then F_max / I_max elseif force_abs < F_N then F_N / I_N else F_N / I_N - Delta_k_F_Force * (force_abs - F_N);
        k_F_0_20 * c_Temp = if T > 393.15 then F_N / I_N elseif T < 293.15 then k_F_0_20 else F_N / I_N + Delta_k_F_Temp * (393.15 - T);
        k_F = k_F_0_20 * c_Temp * c_Force;
        connect(internalSupport.flange, support) annotation(Line(points = {{-80, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(internalSupport.flange, fixed.flange) annotation(Line(points = {{-80, 0}, {-80, -10}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(defaultComponentName = "emf", Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Line(points=  {{0, 90}, {0, 40}}, color=  {0, 0, 255}), Ellipse(extent=  {{-40, 40}, {40, -40}}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, lineColor=  {0, 0, 255}), Line(points=  {{0, -90}, {0, -40}}, color=  {0, 0, 255}), Text(extent=  {{0, -50}, {199, -90}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(visible=  not useSupport, points=  {{-90, -40}, {-30, -40}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-90, -60}, {-70, -40}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-70, -60}, {-50, -40}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-50, -60}, {-30, -40}}, color=  {0, 0, 0}), Line(visible=  not useSupport, points=  {{-60, -40}, {-60, -20}}, color=  {0, 0, 0}), Line(points=  {{-70, -20}, {70, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None, thickness=  0.5), Line(visible=  not useSupport, points=  {{50, 0}, {40, 0}}, color=  {0, 0, 0})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points=  {{-17, 95}, {-20, 85}, {-23, 95}, {-17, 95}}, lineColor=  {160, 160, 164}, fillColor=  {160, 160, 164},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-20, 110}, {-20, 85}}, color=  {160, 160, 164}), Text(extent=  {{-40, 110}, {-30, 90}}, lineColor=  {160, 160, 164}, textString=  "i"), Line(points=  {{9, 75}, {19, 75}}, color=  {192, 192, 192}), Line(points=  {{-20, -110}, {-20, -85}}, color=  {160, 160, 164}), Polygon(points=  {{-17, -100}, {-20, -110}, {-23, -100}, {-17, -100}}, lineColor=  {160, 160, 164}, fillColor=  {160, 160, 164},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-40, -110}, {-30, -90}}, lineColor=  {160, 160, 164}, textString=  "i"), Line(points=  {{8, -79}, {18, -79}}, color=  {192, 192, 192}), Line(points=  {{14, 80}, {14, 70}}, color=  {192, 192, 192})}), Documentation(info = "<html>
<p>EMF transforms electrical energy into rotational mechanical energy. It is used as basic building block of an electrical motor. The mechanical connector flange can be connected to elements of the Modelica.Mechanics.Rotational library. flange.tau is the cut-torque, flange.phi is the angle at the rotational connection.</p>
</html>", revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Martin Otter<br> initially implemented<br>
       </li>
</ul>
</html>"));
      end EMF_linear;

      model MotorResistor "Ideal linear electrical resistor for the d-q system"
        parameter Modelica.SIunits.Resistance R_ref(start = 1)
          "Resistance at temperature T_ref";
        parameter Modelica.SIunits.Temperature T_ref = 293.15
          "Reference temperature";
        parameter Modelica.SIunits.LinearTemperatureCoefficient alpha = 0
          "Temperature coefficient of resistance (R_actual = R*(1 + alpha*(T_heatPort - T_ref))";
        Modelica.SIunits.Resistance R_actual
          "Actual resistance = R*(1 + alpha*(T_heatPort - T_ref))";
        Interfaces.Pin p annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        SI.Current i_d(start = 0) "Current d-coordinates";
        SI.Current i_q(start = 0) "Current q-coordinates";
        SI.Voltage u_d(start = 0) "Voltage d-coordinates";
        SI.Voltage u_q(start = 0) "Voltage q-coordinates";
        SI.Temperature T(start = T_ref) "Temperature of resistance";
        SI.Power LossPower(start = 0) "Lost power at winding resistances";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(Placement(transformation(extent = {{-10, 30}, {10, 50}}), iconTransformation(extent = {{-10, 30}, {10, 50}})));
      equation
        //assert((1 + alpha*(T - T_ref)) >= Modelica.Constants.eps, "Temperature outside scope of model!");
        R_actual = min(R_ref * (1 + (T - T_ref) * alpha), 3 * R_ref);
        u_q = p.u_q - n.u_q;
        u_d = p.u_d - n.u_d;
        i_d = p.i_d;
        i_q = p.i_q;
        0 = p.i_q + n.i_q;
        0 = p.i_d + n.i_d;
        0 = p.w - n.w;
        u_q = R_actual * i_q;
        u_d = R_actual * i_d;
        LossPower = 3 * R_actual * (i_d ^ 2 + i_q ^ 2);
        heatPort.Q_flow = -LossPower;
        heatPort.T = T;
        annotation(Documentation(info = "<html>
<p>The linear resistor connects the branch voltage <i>v</i> with the branch current <i>i</i> by <i>i*R = v</i>. The Resistance <i>R</i> is allowed to be positive, zero, or negative.</p>
</html>", revisions = "<html>
<ul>
<li><i> August 07, 2009   </i>
       by Anton Haumer<br> temperature dependency of resistance added<br>
       </li>
<li><i> March 11, 2009   </i>
       by Christoph Clauss<br> conditional heat port added<br>
       </li>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Rectangle(extent=  {{-70, 30}, {70, -30}}, lineColor=  {0, 0, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-90, 0}, {-70, 0}}, color=  {0, 0, 255}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 255}), Text(extent=  {{-152, -29}, {148, -69}}, textString=  "%name", lineColor=  {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Rectangle(extent=  {{-70, 30}, {70, -30}}, lineColor=  {0, 0, 255}), Line(points=  {{-96, 0}, {-70, 0}}, color=  {0, 0, 255}), Line(points=  {{70, 0}, {96, 0}}, color=  {0, 0, 255})}));
      end MotorResistor;

      model MotorInductance
        "Ideal linear electrical inductance for the d-q system"

        Interfaces.Pin p annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        parameter SI.Inductance L_D "Rotating field inductance";
        SI.AngularVelocity w "Angular velocity of current";
        SI.Current i_d(start = 0) "Current d-coordinates";
        SI.Current i_q(start = 0) "Current q-coordinates";
        SI.Voltage u_d(start = 0) "Voltage d-coordinates";
        SI.Voltage u_q(start = 0) "Voltage q-coordinates";
      equation
        u_d = -w * L_D * i_q;
        //+ L_D*der(i_d);
        u_q = w * L_D * i_d + L_D * der(i_q);
        u_q = p.u_q - n.u_q;
        u_d = p.u_d - n.u_d;
        i_q = p.i_q;
        i_d = p.i_d;
        w = p.w;
        0 = p.i_q + n.i_q;
        0 = p.i_d + n.i_d;
        0 = p.w - n.w;
        annotation(Documentation(info = "<html>
<p>The linear inductor connects the branch voltage <i>v</i> with the branch current <i>i</i> by <i>v = L * di/dt</i>. The Inductance <i>L</i> is allowed to be positive, zero, or negative.</p>
</html>", revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Ellipse(extent=  {{-60, -15}, {-30, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{-30, -15}, {0, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{0, -15}, {30, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{30, -15}, {60, 15}}, lineColor=  {0, 0, 255}), Rectangle(extent=  {{-60, -30}, {60, 0}}, lineColor=  {255, 255, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{60, 0}, {90, 0}}, color=  {0, 0, 255}), Line(points=  {{-90, 0}, {-60, 0}}, color=  {0, 0, 255}), Text(extent=  {{-152, 79}, {148, 39}}, textString=  "%name", lineColor=  {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Ellipse(extent=  {{-60, -15}, {-30, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{-30, -15}, {0, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{0, -15}, {30, 15}}, lineColor=  {0, 0, 255}), Ellipse(extent=  {{30, -15}, {60, 15}}, lineColor=  {0, 0, 255}), Rectangle(extent=  {{-60, -30}, {60, 0}}, lineColor=  {255, 255, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{60, 0}, {96, 0}}, color=  {0, 0, 255}), Line(points=  {{-96, 0}, {-60, 0}}, color=  {0, 0, 255})}));
      end MotorInductance;

      model MotorGround "Ground node for d-q system"

        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{-10, 90}, {10, 110}}), iconTransformation(extent = {{-10, 90}, {10, 110}})));
      equation
        p.u_d = 0;
        p.u_q = 0;
        annotation(Documentation(info = "<html>
<p>Ground of an electrical circuit. The potential at the ground node is zero. Every electrical circuit has to contain at least one ground object.</p>
</html>", revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Line(points=  {{-60, 50}, {60, 50}}, color=  {0, 0, 255}), Line(points=  {{-40, 30}, {40, 30}}, color=  {0, 0, 255}), Line(points=  {{-20, 10}, {20, 10}}, color=  {0, 0, 255}), Line(points=  {{0, 90}, {0, 50}}, color=  {0, 0, 255}), Text(extent=  {{-144, -19}, {156, -59}}, textString=  "%name", lineColor=  {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Line(points=  {{-60, 50}, {60, 50}}, thickness=  0.5, color=  {0, 0, 255}), Line(points=  {{-40, 30}, {40, 30}}, thickness=  0.5, color=  {0, 0, 255}), Line(points=  {{-20, 10}, {20, 10}}, thickness=  0.5, color=  {0, 0, 255}), Line(points=  {{0, 96}, {0, 50}}, thickness=  0.5, color=  {0, 0, 255}), Text(extent=  {{-24, -38}, {22, -6}}, textString=  "p.v=0", lineColor=  {0, 0, 255})}));
      end MotorGround;

      model Friction
        "Model of angular velocity dependent friction losses according to Gross, Hamann and Wiegärtner (2006)"
        extends Modelica.Electrical.Machines.Interfaces.FlangeSupport;
        constant Real eps = 0.1;
        parameter SI.AngularVelocity w_ref = 10
          "Angular velocity after which the friction torque can is approx. equal to friction value";
        parameter Real C_F
          "Constant relating friction torque and square root of angular velocity";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
          "Heat port to model heat flow"                                                            annotation(Placement(transformation(origin = {-100, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 270), iconTransformation(extent = {{10, -10}, {-10, 10}}, rotation = 270, origin = {-100, 0})));
      equation
        tau = -smooth(1, if abs(w) <= eps then 0 else C_F * abs(w) ^ 0.5);
        //tau = -2 / Constants.pi * C_F * abs(w)^0.5 * atan(w / w_ref);
        heatPort.Q_flow = tau * w;
        annotation(Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent=  {{-60, 60}, {60, -60}}, lineColor=  {0, 0, 0}, fillColor=  {175, 175, 175},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-50, 50}, {50, -50}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-12, 50}, {8, 30}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{-10, -30}, {10, -50}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{24, -10}, {44, -30}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{22, 34}, {42, 14}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{-44, 30}, {-24, 10}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{-44, -12}, {-24, -32}}, lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Sphere, fillColor=  {135, 135, 135}), Ellipse(extent=  {{-30, 30}, {30, -30}}, lineColor=  {0, 0, 0}, fillColor=  {175, 175, 175},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-20, 20}, {20, -20}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-260, 101}, {40, 61}}, lineColor=  {0, 0, 255}, textString=  "friction")}), Documentation(info = "<html>
<p>
The friction losses are considered by the equations
</p>
<pre>
   tau / tauRef = (+w / wRef) ^ power_w    for w &gt; +wLinear
 - tau / tauRef = (-w / wRef) ^ power_w    for w &lt; -wLinear
</pre>
<p>
with
</p>
<pre>
  tauRef * wRef = PRef
</pre>
<p>
being the friction torque at the referenc angular velocity
<code>wRef</code>. The exponent <code>power_w</code> is
approximately 1.5 for axial ventilation and approximately 2.0 for radial ventilation.
</p>
<p>
For stability reasons the friction torque <code>tau</code> is approximated by a linear curve
</p>
<pre>
  tau / tauLinear = w / wLinear
</pre>
<p>
with
</p>
<pre>
  tauLinear = tauRef*(wLinear/wRef) ^ power_w
</pre>
<p>
in the range <code> -wLinear &le; w &le; wLinear</code> with <code>wLinear = 0.001 * wRef</code>. The relationship of torque
and angular velocity is depicted in Fig. 1
</p>
<table border=0 cellspacing=0 cellpadding=1>
  <tr><td> <img src=\"modelica://Modelica/Resources/Images/Electrical/Machines/frictiontorque.png\"> </td>
  </tr>
  <tr><td> <b> Fig. 1: </b>Friction loss torque versus angular velocity for <code>power_w = 2</code></td>
  </tr>
</table>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Electrical.Machines.Losses.FrictionParameters\">FrictionParameters</a>
</p>
<p>
If it is desired to neglect friction losses, set <code>frictionParameters.PRef = 0</code> (this is the default).
</p>
</html>"));
      end Friction;

      model MotorCurrentSensor "Sensor to measure the current in a branch"
        extends Modelica.Icons.RotationalSensor;
        Modelica.Blocks.Interfaces.RealOutput i_actual
          "current in the branch from p to n as output signal"                                              annotation(Placement(transformation(origin = {0, -100}, extent = {{10, -10}, {-10, 10}}, rotation = 90), iconTransformation(extent = {{10, -10}, {-10, 10}}, rotation = 90, origin = {20, -99})));
        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{-110, -11}, {-90, 9}}), iconTransformation(extent = {{-110, -11}, {-90, 9}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Interfaces.RealOutput i_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{87, -50}, {107, -30}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {40, -99})));
        Modelica.Blocks.Interfaces.RealOutput i_rms
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{88, -80}, {108, -60}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {0, -99})));
        Modelica.Blocks.Interfaces.RealOutput i_q
          "current in the branch from p to n as output signal"                                         annotation(Placement(transformation(origin = {-20, -99}, extent = {{10, -10}, {-10, 10}}, rotation = 90), iconTransformation(extent = {{10, -10}, {-10, 10}}, rotation = 90, origin = {-20, -99})));
        Modelica.Blocks.Interfaces.RealOutput i_d
          "current in the branch from p to n as output signal"                                         annotation(Placement(transformation(origin = {-40, -99}, extent = {{10, -10}, {-10, 10}}, rotation = 90), iconTransformation(extent = {{10, -10}, {-10, 10}}, rotation = 90, origin = {-40, -99})));
      protected
        Modelica.Blocks.Interfaces.RealInput u1
          "Connector of Real input signal"                                       annotation(Placement(transformation(extent = {{5, -79}, {23, -61}})));
      public
        HelpBlocks.RootMeanSquareValue rootMeanSquareValue annotation(Placement(transformation(extent = {{42, -80}, {62, -60}})));
        HelpBlocks.Max max
          annotation (Placement(transformation(extent={{41,-50},{61,-30}})));
      equation
        p.u_q = n.u_q;
        p.u_d = n.u_d;
        p.w = n.w;
        0 = p.i_q + n.i_q;
        0 = p.i_d + n.i_d;
        i_actual = sqrt(p.i_q ^ 2 + p.i_d ^ 2);
        i_q = p.i_q;
        i_d = p.i_d;
        connect(i_rms, i_rms) annotation(Line(points = {{98, -70}, {98, -94.5}, {98, -94.5}, {98, -70}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(i_actual, u1) annotation(Line(points = {{0, -100}, {10, -100}, {10, -70}, {14, -70}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(u1, rootMeanSquareValue.u) annotation(Line(points = {{14, -70}, {40, -70}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue.y, i_rms) annotation(Line(points = {{63, -70}, {98, -70}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(max.y, i_max) annotation (Line(
            points={{62,-40},{97,-40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(max.u, rootMeanSquareValue.u) annotation (Line(
            points={{39,-40},{28,-40},{28,-70},{40,-70}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Text(extent=  {{-29, -11}, {30, -70}}, lineColor=  {0, 0, 0}, textString=  "A"), Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Text(extent=  {{-150, 80}, {150, 120}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Line(points=  {{0, -90}, {0, -70}}, color=  {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}},                                                                                                    grid = {1, 1}), graphics={  Text(extent=  {{-153, 79}, {147, 119}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{-70, 0}, {-96, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {96, 0}}, color=  {0, 0, 0}), Line(points=  {{0, -90}, {0, -70}}, color=  {0, 0, 255})}), Documentation(revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>", info = "<html>
<p>The current  sensor converts the current flowing between the two connectors into a real valued signal. The two connectors are in the sensor connected like a short cut. The sensor has to be placed within an electrical connection in series.  It does not influence the current sum at the connected nodes. Therefore, the electrical behavior is not influenced by the sensor.</p>
</html>"));
      end MotorCurrentSensor;

      model MotorVoltageSensor "Sensor to measure the voltage between two pins"
        extends Modelica.Icons.RotationalSensor;
        Modelica.Blocks.Interfaces.RealOutput u_link_actual(start = 0)
          "Voltage between pin p and n (= p.v - n.v) as output signal"                                                              annotation(Placement(transformation(origin={10,-100},   extent = {{10, -10}, {-10, 10}}, rotation = 90), iconTransformation(extent = {{10, -10}, {-10, 10}}, rotation = 90, origin={10,-100})));
        SI.Voltage u_d(start = 0);
        SI.Voltage u_q(start = 0);
        Interfaces.PositivePin p annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Interfaces.RealOutput u_link_max
          "Connector of Real output signal"                                                annotation(Placement(transformation(extent = {{62, -69}, {82, -49}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin={30,-100})));
      protected
        Modelica.Blocks.Interfaces.RealInput u1
          "Connector of Real input signal"                                       annotation(Placement(transformation(extent = {{0, -69}, {21, -49}})));
      public
        Modelica.Blocks.Interfaces.RealOutput u_q_out
          "Connector of Real output signal" annotation (Placement(transformation(
                extent={{62,-69},{82,-49}}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,-100})));
        Modelica.Blocks.Interfaces.RealOutput u_d_out(start=0)
          "Voltage between pin p and n (= p.v - n.v) as output signal" annotation (
            Placement(transformation(
              origin={10,-100},
              extent={{10,-10},{-10,10}},
              rotation=90), iconTransformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-30,-100})));
        HelpBlocks.Max max
          annotation (Placement(transformation(extent={{29,-69},{49,-49}})));
      equation
        p.i_q = 0;
        n.i_q = 0;
        p.i_d = 0;
        n.i_d = 0;
        u_d = p.u_d - n.u_d;
        u_q = p.u_q - n.u_q;
        u_link_actual = sqrt(3) * sqrt(u_d ^ 2 + u_q ^ 2);
        u_d_out = u_d;
        u_q_out = u_q;
        connect(u_link_actual, u1) annotation(Line(points={{10,-100},{8,-100},{8,-59},
                {10.5,-59}},                                                                               color = {0, 0, 127}, smooth = Smooth.None));
        connect(u1, max.u) annotation (Line(
            points={{10.5,-59},{27,-59}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(max.y, u_link_max) annotation (Line(
            points={{50,-59},{72,-59}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Icon(coordinateSystem(preserveAspectRatio=false,   extent={{-100,-100},
                  {100,100}},                                                                              grid = {1, 1}), graphics={  Text(extent=  {{-29, -11}, {30, -70}}, lineColor=  {0, 0, 0}, textString=  "V"), Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Line(points=  {{0, -90}, {0, -70}}, color=  {0, 0, 255}), Text(extent=  {{-150, 80}, {150, 120}}, textString=  "%name", lineColor=  {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}},                                                                                                    grid = {1, 1}), graphics={  Line(points=  {{-70, 0}, {-96, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {96, 0}}, color=  {0, 0, 0}), Line(points=  {{0, -90}, {0, -70}}, color=  {0, 0, 255})}), Documentation(revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>", info = "<html>
<p>The voltage  sensor converts the voltage between the two connectors into a real valued signal. It does not influence the current sum at the nodes in between the voltage is measured, therefore, the electrical behavior is not influenced by the sensor.</p>
</html>"));
      end MotorVoltageSensor;

      model PowerSensor "Sensor to measure the power"
        Modelica.Blocks.Interfaces.RealOutput power annotation(Placement(transformation(origin = {-80, -110}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
        Modelica.Blocks.Interfaces.RealOutput energy annotation(Placement(transformation(origin = {-80, -110}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
        Modelica.Blocks.Interfaces.RealOutput power_average annotation(Placement(transformation(origin = {-80, -110}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
        Interfaces.Pin p_c annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.NegativePin n_c annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Interfaces.PositivePin p_v annotation(Placement(transformation(extent = {{-10, 86}, {10, 106}}), iconTransformation(extent = {{-10, 86}, {10, 106}})));
        Interfaces.NegativePin n_v annotation(Placement(transformation(extent = {{-10, -110}, {10, -90}}), iconTransformation(extent = {{-10, -110}, {10, -90}})));
      protected
        parameter Real startTime(fixed = false);
      initial equation
        startTime = time;
        energy = 0;
      equation
        power = 3 * p_c.i_q * (p_v.u_q - n_v.u_q) + 3 * p_c.i_d * (p_v.u_d - n_v.u_d);
        der(energy) = power;
        if time <= startTime then
          power_average = energy;
        else
          power_average = energy / (time - startTime);
        end if;
        0 = p_c.u_q - n_c.u_q;
        0 = p_c.u_d - n_c.u_d;
        0 = p_c.w - n_c.w;
        0 = p_c.i_q + n_c.i_q;
        0 = p_c.i_d + n_c.i_d;
        p_v.i_q = 0;
        n_v.i_q = 0;
        p_v.i_d = 0;
        n_v.i_d = 0;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Ellipse(extent=  {{-70, 70}, {70, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 100}, {0, 70}}, color=  {0, 0, 255}), Line(points=  {{0, -70}, {0, -100}}, color=  {0, 0, 255}), Line(points=  {{-80, -100}, {-80, 0}}, color=  {0, 0, 255}), Line(points=  {{-100, 0}, {100, 0}}, color=  {0, 0, 255}), Text(extent=  {{150, 120}, {-150, 160}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{0, 70}, {0, 40}}, color=  {0, 0, 0}), Line(points=  {{22.9, 32.8}, {40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-22.9, 32.8}, {-40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{37.6, 13.7}, {65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-37.6, 13.7}, {-65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{0, 0}, {9.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-0.48, 31.6}, {18, 26}, {18, 57.2}, {-0.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-5, 5}, {5, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-29, -11}, {30, -70}}, lineColor=  {0, 0, 0}, textString=  "P")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics), Documentation(info = "<html>
<p>This power sensor measures instantaneous electrical power of a singlephase system and has a separated voltage and current path. The pins of the voltage path are pv and nv, the pins of the current path are pc and nc. The internal resistance of the current path is zero, the internal resistance of the voltage path is infinite.</p>
</html>", revisions = "<html>
<ul>
<li><i>January 12, 2006</i> by Anton Haumer implemented</li>
</ul>
</html>"));
      end PowerSensor;

      model v_p_controller "combined speed and position controller"
        parameter Real K_V "(m/min)/mm";
        parameter Real K_pn "A/(m/s)";
        parameter Real T_nn "s";
        parameter Real delay(displayUnit = "ms");
        Modelica.Blocks.Continuous.Integrator integrator annotation(Placement(transformation(extent = {{-60, 70}, {-80, 90}})));
        Modelica.Blocks.Continuous.PI PI_v(k = K_pn, T = T_nn) annotation(Placement(transformation(extent = {{-6, -10}, {14, 10}})));
        HelpBlocks.Feedback_mirror v_feedback annotation(Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 180, origin = {-28, 0})));
        Modelica.Blocks.Math.Gain gain(k = K_V * 16.67) annotation(Placement(transformation(extent = {{-68, -10}, {-48, 10}})));
        HelpBlocks.Feedback_mirror s_feedback annotation(Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 180, origin = {-88, 0})));
        Modelica.Blocks.Sources.Step step(height = 100, startTime = 0.01) annotation(Placement(transformation(extent = {{200, 0}, {180, 20}})));
        Modelica.Blocks.Nonlinear.FixedDelay fixedDelay(delayTime = delay) annotation(Placement(transformation(extent = {{-4, 70}, {-24, 90}})));
        Modelica.Blocks.Interfaces.RealInput x_s annotation(Placement(transformation(extent = {{-140, -20}, {-100, 20}}), iconTransformation(extent = {{-140, -20}, {-100, 20}})));
        Modelica.Blocks.Interfaces.RealInput v_a
          "Connector of Real input signal"                                        annotation(Placement(transformation(extent = {{40, 74}, {0, 114}}), iconTransformation(extent = {{20, -20}, {-20, 20}}, rotation = 90, origin = {0, 120})));
        Modelica.Blocks.Interfaces.RealOutput i_s
          "Connector of Real output signal"                                         annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
      equation
        connect(PI_v.u, v_feedback.y) annotation(Line(points = {{-8, 0}, {-19, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(gain.y, v_feedback.u1) annotation(Line(points = {{-47, 0}, {-36, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(s_feedback.y, gain.u) annotation(Line(points = {{-79, 0}, {-70, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(integrator.y, s_feedback.u2) annotation(Line(points = {{-81, 80}, {-88, 80}, {-88, 8}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(fixedDelay.y, integrator.u) annotation(Line(points = {{-25, 80}, {-58, 80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(fixedDelay.y, v_feedback.u2) annotation(Line(points = {{-25, 80}, {-28, 80}, {-28, 8}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(s_feedback.u1, x_s) annotation(Line(points = {{-96, 0}, {-120, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(fixedDelay.u, v_a) annotation(Line(points = {{-2, 80}, {54, 80}, {54, 94}, {20, 94}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(PI_v.y, i_s) annotation(Line(points = {{15, 0}, {110, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}), Text(extent = {{-80, 20}, {80, -20}}, lineColor = {0, 0, 0}, textString = "vel & pos 
control"), Text(extent = {{-260, -118}, {260, -160}}, lineColor = {0, 0, 0}, textString = "K_V = %K_V m/min/mm"), Text(extent = {{-260, -162}, {260, -204}}, lineColor = {0, 0, 0}, textString = "K_Pn = %K_pn A/(m/s)"), Text(extent = {{-260, -210}, {260, -252}}, lineColor = {0, 0, 0}, textString = "T_nn = %T_nn s")}));
      end v_p_controller;

      block voltageLimit
        "Limit the range of a voltage signal with variable limits"
        extends Modelica.Blocks.Interfaces.SISO;
        parameter Boolean strict=false
          "= true, if strict limits with noEvent(..)"
          annotation (Evaluate=true, choices(checkBox=true));
        parameter Boolean limitsAtInit=true
          "= false, if limits are ignored during initialization (i.e., y=u)"
          annotation (Evaluate=true, choices(checkBox=true));
        parameter SI.Voltage U_eff_max "Voltage limit of motor module" annotation(Dialog(group = "Limiting data"));
        Modelica.Blocks.Interfaces.RealInput u_d
          "Connector of Real input signal used as maximum of input u" annotation (
            Placement(transformation(extent={{-140,60},{-100,100}}, rotation=0)));

      protected
        Real uMax;
        Real uMin;
        Real limit1;
        Real limit2;

      equation
        if strict then
          uMax = noEvent(max(limit1, limit2));
          uMin = noEvent(min(limit1, limit2));
        else
          uMax = max(limit1, limit2);
          uMin = min(limit1, limit2);
        end if;
        limit2 = - limit1;
        limit1 =  sqrt(max(0,U_eff_max^2/3 - u_d^2));
        if initial() and not limitsAtInit then
           y = u;
           assert(u >= uMin - 0.01*abs(uMin) and
                  u <= uMax + 0.01*abs(uMax),
                 "VariableLimiter: During initialization the limits have been ignored.\n"+
                 "However, the result is that the input u is not within the required limits:\n"+
                 "  u = " + String(u) + ", uMin = " + String(uMin) + ", uMax = " + String(uMax));
        elseif strict then
          y = smooth(0, noEvent(if u > uMax then uMax else if u < uMin then uMin else u));
        else
           y = smooth(0,if u > uMax then uMax else if u < uMin then uMin else u);
        end if;
        annotation (
          Documentation(info="<html>
<p>
The Limiter block passes its input signal as output signal
as long as the input is within the upper and lower
limits specified by the two additional inputs limit1 and
limit2. If this is not the case, the corresponding limit
is passed as output.
</p>
</html>"),       Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={
              Line(points={{0,-90},{0,68}}, color={192,192,192}),
              Line(points={{-90,0},{68,0}}, color={192,192,192}),
              Polygon(
                points={{90,0},{68,-8},{68,8},{90,0}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,-70},{-50,-70},{50,70},{80,70}}, color={0,0,0}),
              Text(
                extent={{-150,150},{150,110}},
                textString="%name",
                lineColor={0,0,255}),
              Line(points={{-100,80},{66,80},{66,70}}, color={0,0,127}),
              Line(points={{-100,-80},{-64,-80},{-64,-70}}, color={0,0,127}),
              Polygon(points={{-64,-70},{-66,-74},{-62,-74},{-64,-70}}, lineColor={
                    0,0,127}),
              Polygon(points={{66,70},{64,74},{68,74},{66,70}}, lineColor={0,0,127}),
              Polygon(
                points={{0,90},{-8,68},{8,68},{0,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(
                visible=strict,
                points={{50,70},{80,70}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                visible=strict,
                points={{-80,-70},{-50,-70}},
                color={255,0,0},
                smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                  100}}),     graphics={
              Line(points={{0,-60},{0,50}}, color={192,192,192}),
              Polygon(
                points={{0,60},{-5,50},{5,50},{0,60}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-60,0},{50,0}}, color={192,192,192}),
              Polygon(
                points={{60,0},{50,-5},{50,5},{60,0}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-50,-40},{-30,-40},{30,40},{50,40}}, color={0,0,0}),
              Text(
                extent={{46,-6},{68,-18}},
                lineColor={128,128,128},
                textString="inPort"),
              Text(
                extent={{-30,70},{-5,50}},
                lineColor={128,128,128},
                textString="outPort"),
              Text(
                extent={{-66,-40},{-26,-20}},
                lineColor={128,128,128},
                textString="uMin"),
              Text(
                extent={{30,20},{70,40}},
                lineColor={128,128,128},
                textString="uMax"),
              Line(points={{-100,80},{40,80},{40,40}}, color={0,0,127}),
              Polygon(points={{40,40},{35,50},{45,50},{40,40}}, lineColor={0,0,127}),
              Polygon(points={{-40,-40},{-45,-50},{-35,-50},{-40,-40}}, lineColor={
                    0,0,127})}));
      end voltageLimit;
    end Components;

    package Interfaces
      connector Pin "connector for motor d-q system"
        Modelica.SIunits.Voltage u_d(start = 0) "Potential at the pin";
        Modelica.SIunits.Voltage u_q(start = 0) "Potential at the pin" annotation(unassignedMessage = "An electrical potential cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        flow Modelica.SIunits.Current i_d "Current flowing into the pin";
        flow Modelica.SIunits.Current i_q "Current flowing into the pin" annotation(unassignedMessage = "An electrical current cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        Modelica.SIunits.AngularVelocity w "Angular speed of current";
        annotation(defaultComponentName = "pin", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-40, 40}, {40, -40}}, lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-160, 110}, {40, 50}}, lineColor=  {0, 0, 255}, textString=  "%name")}), Documentation(revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>", info = "<html>
<p>Pin is the basic electric connector. It includes the voltage which consists between the pin and the ground node. The ground node is the node of (any) ground device (Modelica.Electrical.Basic.Ground). Furthermore, the pin includes the current, which is considered to be <b>positive</b> if it is flowing at the pin<b> into the device</b>.</p>
</html>"));
      end Pin;

      connector PositivePin "Positive pin of an electric component"
        Modelica.SIunits.Voltage u_d(start = 0) "Potential at the pin";
        Modelica.SIunits.Voltage u_q(start = 0) "Potential at the pin" annotation(unassignedMessage = "An electrical potential cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        flow Modelica.SIunits.Current i_d "Current flowing into the pin";
        flow Modelica.SIunits.Current i_q "Current flowing into the pin" annotation(unassignedMessage = "An electrical current cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        Modelica.SIunits.AngularVelocity w "Potential at the pin";
        annotation(defaultComponentName = "pin_p", Documentation(info = "<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>", revisions = "<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-40, 40}, {40, -40}}, lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-160, 110}, {40, 50}}, lineColor=  {0, 0, 255}, textString=  "%name")}));
      end PositivePin;

      connector NegativePin "Negative pin of an electric component"
        Modelica.SIunits.Voltage u_d(start = 0) "Potential at the pin";
        Modelica.SIunits.Voltage u_q(start = 0) "Potential at the pin" annotation(unassignedMessage = "An electrical potential cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        flow Modelica.SIunits.Current i_d "Current flowing into the pin";
        flow Modelica.SIunits.Current i_q "Current flowing into the pin" annotation(unassignedMessage = "An electrical current cannot be uniquely calculated.
                                                                        The reason could be that
                                                                        - a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
                                                                          to define the zero potential of the electrical circuit, or
                                                                        - a connector of an electrical component is not connected.");
        Modelica.SIunits.AngularVelocity w "Potential at the pin";
        annotation(defaultComponentName = "pin_n", Documentation(info = "<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>", revisions = "<html>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-40, 40}, {40, -40}}, lineColor=  {0, 0, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-40, 110}, {160, 50}}, textString=  "%name", lineColor=  {0, 0, 255})}));
      end NegativePin;
    end Interfaces;
  end FeedDriveMotor;

  package LinearActuators
    package BallScrew
      type BearingConfiguration = enumeration(
          fixed_free "Bearing factor = 34.9 m/min",
          supported_supported "Bearing factor = 108 m/min",
          fixed_supported "Bearing factor = 168 m/min",
          fixed_fixed "Bearing factor = 242 m/min")
        "Enumeration defining spindle support bearing configuration"                                                                                                     annotation(Evaluate = true);

      model BallScrewDrive
        "Model of a ball screw drive according to Gross, Hamann and Wiegärtner (2006). Can be parameterized by typical supplier data such as from A. Mannesmann"
        parameter SI.Diameter d(displayUnit = "mm")
          "Nominal diameter of spindle"                                           annotation(Dialog(group = "Geometric data"));
        parameter SI.Length P(displayUnit = "mm") "Pitch of spindle" annotation(Dialog(group = "Geometric data"));
        parameter SI.Length l_S = 1 "Length of spindle" annotation(Dialog(group = "Geometric data"));
        parameter SI.Length l_S_u = 0.9 * l_S "Unsupported length of spindle" annotation(Dialog(group = "Geometric data", tab = "Calculated"));
        parameter SI.Density density = 7900 "Density of spindle material" annotation(Dialog(group = "Physical data"));
        parameter Modelica.SIunits.Mass m = (d / 2) ^ 2 * Constants.pi * l_S * density
          "Mass of Spindle"                                                                              annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter SI.Conversions.NonSIunits.Angle_deg phi = SI.Conversions.to_deg(atan(P / (Constants.pi * d)))
          "Angle of lead"                                                                                                     annotation(Dialog(group = "Geometric data", tab = "Calculated"));
        parameter SI.Conversions.NonSIunits.Angle_deg rho = 0.2
          "Angle of friction"                                                       annotation(Dialog(group = "Geometric data"));
        parameter Real eta = tan(SI.Conversions.from_deg(phi)) / tan(SI.Conversions.from_deg(phi + rho))
          "Value of efficiency of nut"                                                                                                annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter SI.MomentOfInertia J = 0.5 * m * d ^ 2 / 4
          "Moment of inertia of the spindle"                                                    annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter SI.Force C_am = 50000 "Dynamic load rating of nut" annotation(Dialog(group = "Engineering data"));
        parameter SI.Force C_0am = 80000 "Static load rating of nut" annotation(Dialog(group = "Engineering data"));
        parameter Real preload_factor = 0.07
          "Preload of factor for nut. preload_factor = F_pr/C_am"                                    annotation(Dialog(group = "Engineering data"));
        parameter SI.TranslationalSpringConstant R_nu = 0.5 * 10 ^ 9
          "Nut rigidity on nut flange"                                                            annotation(Dialog(group = "Physical data"));
        parameter SI.TranslationalSpringConstant R_S = 0.09 * 10 ^ 9
          "Spindle rod rigidity per m"                                                            annotation(Dialog(group = "Physical data"));
        parameter SI.TranslationalSpringConstant R_ax = 1 / (1 / R_nu + l_S_u / one_meter / (R_S * k))
          "Axial rigidity"                                                                                              annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter Real mu = 0.006 "Friction coefficient of bearing" annotation(Dialog(group = "Physical data"));
        parameter SI.Force F_pr_Bearing = 0 "Preload force of bearing" annotation(Dialog(group = "Engineering data"));
        parameter SI.Diameter d_B = 1.6 * d "Diameter of bearing" annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        parameter BearingConfiguration bearingConfiguration = BearingConfiguration.fixed_supported
          "Bearing method of spindle"                                                                                          annotation(Dialog(group = "Engineering data"), Evaluate = true, choices(choice = LinearActuators.BallScrew.BearingConfiguration.fixed_free
              "fixed-free, k_L = 34.9e3 m/min",                                                                                                    choice = LinearActuators.BallScrew.BearingConfiguration.supported_supported
              "supported-supported, k_L = 108e3 m/min",                                                                                                    choice = LinearActuators.BallScrew.BearingConfiguration.fixed_supported
              "fixed-supported, k_L = 168e3 m/min",                                                                                                    choice = LinearActuators.BallScrew.BearingConfiguration.fixed_fixed
              "fixed-fixed, k_L = 242e3 m/min"));
        parameter SI.Force F_pr = preload_factor * C_am "Preload force of nut" annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        parameter SI.Mass m_ref = 1000 "Reference mass for eigenfrequency" annotation(Dialog(group = "Engineering data"));
        parameter SI.Frequency f_bs = 1 / (2 * Constants.pi) * sqrt(R_ax / m_ref)
          "Eigenfrequncy of ball screw"                                                                         annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        parameter SI.Frequency f_min = 50
          "Minimum eigenfrequncy of ball screw and spindle. Constraint f_bs > f_min"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_krit = k_L * d / l_S_u ^ 2
          "critcial rotational speed of spindle"                                                                                    annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        parameter SI.Conversions.NonSIunits.Time_hour L_h_min = 10000
          "Minimum life expectancy of nut"                                                             annotation(Dialog(group = "Limiting data"));
        parameter Real DN_perm(unit = "mm.min-1") = 160000
          "Permissible DN value of nut (diamater*max_speed)"                                                  annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_b = k_buckling * d ^ 4 / l_S_u ^ 2
          "Buckling force"                                                       annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        constant SI.Length one_meter = 1;
        parameter Real S_n_krit = 1.25 "Safety factor for n_krit" annotation(Dialog(tab = "Safety factors"));
        parameter Real S_F_b = 2 "Safety factor for buckling" annotation(Dialog(tab = "Safety factors"));
        parameter Real S_C_0am = 3 "Safety factor for static load rate" annotation(Dialog(tab = "Safety factors"));
        SI.Frequency con_f "Constraint for eigenfrequncy of axis";
        SI.Conversions.NonSIunits.Time_hour con_L_h
          "Constraint for life expectancy. The value at t_end is relevant when simulation terminates.";
        SI.Force con_F_pr "Constraint for preload";
        SI.Conversions.NonSIunits.AngularVelocity_rpm con_n_krit
          "Constraint for critcal speed";
        Real con_DN(unit = "mm.min-1");
        SI.Force con_F_b "Constraint for buckling force";
        SI.Force con_C_0am "Constraint for static load rating";
        Real util_f;
        Real util_L_h;
        Real util_F_pr;
        Real util_n_krit;
        Real util_DN;
        Real util_F_b;
        Real util_C_0am;
      protected
        parameter Integer n_Bearing_ax = if bearingConfiguration == BearingConfiguration.supported_supported or bearingConfiguration == BearingConfiguration.fixed_fixed then 2 else 1
          "number of bearing with axial load";
        parameter Integer k = if bearingConfiguration == BearingConfiguration.supported_supported or bearingConfiguration == BearingConfiguration.fixed_fixed then 4 else 1;
        parameter Real k_L(unit = "m.min-1") = if bearingConfiguration == BearingConfiguration.supported_supported then 108 * 10 ^ 3 else if bearingConfiguration == BearingConfiguration.fixed_supported then 168 * 10 ^ 3 else if bearingConfiguration == BearingConfiguration.fixed_fixed then 242 * 10 ^ 3 else 34.9 * 10 ^ 3;
        parameter Real k_buckling(unit = "N/m2") = if bearingConfiguration == BearingConfiguration.supported_supported then 5.6 * 10 ^ 10 else if bearingConfiguration == BearingConfiguration.fixed_supported then 11.2 * 10 ^ 10 else if bearingConfiguration == BearingConfiguration.fixed_fixed then 22.4 * 10 ^ 10 else 1.4 * 10 ^ 10;
      public
        BallScrewNut ballScrewNut(eta = eta, P = P, C_am = C_am, F_pr = F_pr) annotation(Placement(transformation(extent = {{40, -10}, {60, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flangeT
          "(right) driven flange (flange axis directed out of cut plane)"                                                            annotation(Placement(transformation(extent = {{92, -10}, {112, 10}}), iconTransformation(extent = {{92, -10}, {112, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flangeR
          "Left flange of shaft"                                                         annotation(Placement(transformation(extent = {{-100, -10}, {-80, 10}}), iconTransformation(extent = {{-100, -10}, {-80, 10}})));
        Modelica.Mechanics.Rotational.Components.Inertia spindle(J = J) annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        RotationalComponents.BasicComponents.Bearing bearing(mu = mu * n_Bearing_ax, d = d_B, F_pr = F_pr_Bearing) annotation(Placement(transformation(extent = {{-50, -40}, {-30, -20}})));
        Modelica.Mechanics.Rotational.Components.Fixed fixed annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-40, -56})));
        Modelica.Mechanics.Translational.Interfaces.Support supportT1
          "Translational support/housing of component"                                                             annotation(Placement(transformation(extent = {{-110, -40}, {-90, -20}}), iconTransformation(extent = {{-112, -80}, {-92, -60}})));
      initial equation
        util_f = 0;
        util_L_h = 0;
        util_F_pr = 0;
        util_n_krit = 0;
        util_DN = 0;
        util_F_b = 0;
        util_C_0am = 0;
      equation
        con_f = f_min - f_bs;
        con_L_h = L_h_min - ballScrewNut.nominalLifeRating.y;
        con_F_pr = ballScrewNut.force_absMax.y - 2 ^ (3 / 2) * F_pr;
        con_n_krit = S_n_krit * ballScrewNut.speed_absMax.y - n_krit;
        con_DN = ballScrewNut.speed_absMax.y * (d * 1000) - DN_perm;
        con_F_b = S_F_b * ballScrewNut.force_absMax.y - F_b;
        con_C_0am = S_C_0am * ballScrewNut.force_absMax.y - C_0am;
        connect(ballScrewNut.flangeT, flangeT) annotation(Line(points = {{60, 0}, {102, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(spindle.flange_b, ballScrewNut.flangeR) annotation(Line(points = {{0, 0}, {39.6, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(spindle.flange_a, flangeR) annotation(Line(points = {{-20, 0}, {-90, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(bearing.flange, flangeR) annotation(Line(points = {{-40, -20}, {-40, 0}, {-90, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(bearing.support, fixed.flange) annotation(Line(points = {{-40, -40}, {-40, -56}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(bearing.flange_b, ballScrewNut.supportT1) annotation(Line(points = {{-30, -30}, {39.8, -30}, {39.8, -5}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(flangeT, flangeT) annotation(Line(points = {{102, 0}, {102, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(bearing.flange_a, supportT1) annotation(Line(points = {{-50, -30}, {-100, -30}}, color = {0, 127, 0}, smooth = Smooth.None));
        when terminal() then
          util_f = f_min / f_bs;
          util_L_h = L_h_min / ballScrewNut.nominalLifeRating.y;
          util_F_pr = ballScrewNut.force_absMax.y / (2 ^ (3 / 2) * F_pr);
          util_n_krit = S_n_krit * ballScrewNut.speed_absMax.y / n_krit;
          util_DN = ballScrewNut.speed_absMax.y * (d * 1000) / DN_perm;
          util_F_b = S_F_b * ballScrewNut.force_absMax.y / F_b;
          util_C_0am = S_C_0am * ballScrewNut.force_absMax.y / C_0am;
        end when;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{-150, 100}, {150, 60}}, textString = "%name", lineColor = {0, 0, 255}), Rectangle(extent = {{-90, 20}, {90, -20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-50, 20}, {-30, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-30, 20}, {-10, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-10, 20}, {10, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{10, 20}, {30, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{30, 20}, {50, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{50, 20}, {70, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Rectangle(extent = {{-32, 28}, {22, -28}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-12, 40}, {8, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-70, 20}, {-50, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{70, 20}, {90, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{62, 36}, {8, 36}, {102, 36}, {102, 10}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-90, 20}, {-70, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Rectangle(extent = {{-86, 30}, {-66, 20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{-86, -20}, {-66, -30}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{64, 30}, {84, 20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{64, -20}, {84, -30}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Backward), Line(points = {{-76, -30}, {-76, -70}, {-96, -70}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-160, -50}, {160, -86}}, lineColor = {0, 0, 0}, textString = "P = %P")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
      end BallScrewDrive;

      model NominalLifeRating
        "Calculation of nominal life rating for ball screw (see A. Mannesmann)"
        extends Modelica.Blocks.Interfaces.SI2SO;
        SI.Force F_a "Axial operation force";
        SI.Force F_w "Operation force with preload";
        SI.Force F_w_m "Effective operation force";
        parameter SI.Force C_am = 50000 "Dynamic load rating";
        parameter SI.Force F_pr = 0.07 * C_am "Preload force";
        SI.AngularVelocity w "Angular velocity of ball screw";
        SI.AngularVelocity w_m "Average angular velocity";
        Real temp_F(start = 1);
        Real temp_w(start = 1);
        Real L_10 "Number of 10^6 revolution with 90 % life expectancy";
        SI.Conversions.NonSIunits.Time_hour L_h "Life expectancy on 90 % level";
      initial equation
        F_w_m = 0;
        L_10 = 0;
        L_h = 0;
      equation
        F_a = abs(u1);
        w = abs(u2);
        F_w = max(0.5 * F_a + F_pr, F_a);
        der(temp_F) = F_w ^ 3 * w;
        der(temp_w) = w;
        when terminal() then
          F_w_m = (temp_F / temp_w) ^ (1 / 3);
          w_m = temp_w / time;
          L_10 = (C_am / F_w_m) ^ 3;
          L_h = 2 * (L_10 * 10 ^ 6 * 2 * Constants.pi) / (w_m * 3600);
          //assuming load directions with equal distribution
        end when;
        y = L_h;
        annotation(Documentation(info = "
<HTML>
<p>
This blocks computes output <b>y</b> as <i>sum</i> of the
two input signals <b>u1</b> and <b>u2</b>:
</p>
<pre>
    <b>y</b> = k1*<b>u1</b> + k2*<b>u2</b>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   k1= +2, k2= -3

  results in the following equations:

     y = 2 * u1 - 3 * u2
</pre>

</HTML>
"),       Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Text(
                      extent={{-98,-52},{7,-92}},
                      lineColor={0,0,0},
                      textString="%k2"),Text(
                      extent={{-100,90},{5,50}},
                      lineColor={0,0,0},
                      textString="%k1"),Text(
                      extent={{-150,150},{150,110}},
                      textString="%name",
                      lineColor={0,0,255}),Line(points={{-100,60},{-40,60},{-30,
                40}}, color={0,0,255}),Ellipse(extent={{-50,50},{50,-50}},
                lineColor={0,0,255}),Line(points={{-100,-60},{-40,-60},{-30,-40}},
                color={0,0,255}),Line(points={{-15,-25.99},{15,25.99}}, color={
                0,0,0}),Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,127},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Line(points={{50,0},{100,0}},
                color={0,0,255}),Line(points={{-100,60},{-74,24},{-44,24}},
                color={0,0,127}),Line(points={{-100,-60},{-74,-28},{-42,-28}},
                color={0,0,127}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,127}),Line(points={{50,0},{100,0}}, color={0,0,127}),Text(
                      extent={{-38,34},{38,-34}},
                      lineColor={0,0,0},
                      textString="L_h"),Text(
                      extent={{-100,52},{5,92}},
                      lineColor={0,0,0},
                      textString="F_a"),Text(
                      extent={{-100,-52},{5,-92}},
                      lineColor={0,0,0},
                      textString="w")}),
          Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,255},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Text(
                      extent={{-98,-52},{7,-92}},
                      lineColor={0,0,0},
                      textString="%k2"),Text(
                      extent={{-100,90},{5,50}},
                      lineColor={0,0,0},
                      textString="%k1"),Line(points={{-100,60},{-40,60},{-30,40}},
                color={0,0,255}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,255}),Line(points={{-100,-60},{-40,-60},{-30,-40}}, color=
                {0,0,255}),Line(points={{-15,-25.99},{15,25.99}}, color={0,0,0}),
                Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,127},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Line(points={{50,0},{100,0}},
                color={0,0,255}),Line(points={{-100,60},{-74,24},{-44,24}},
                color={0,0,127}),Line(points={{-100,-60},{-74,-28},{-42,-28}},
                color={0,0,127}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,127}),Line(points={{50,0},{100,0}}, color={0,0,127}),Text(
                      extent={{-100,52},{5,92}},
                      lineColor={0,0,0},
                      textString="F_a"),Text(
                      extent={{-100,-52},{5,-92}},
                      lineColor={0,0,0},
                      textString="n"),Text(
                      extent={{-52,-20},{53,20}},
                      lineColor={0,0,0},
                      textString="L_h")}));
      end NominalLifeRating;

      model BallScrewNut
        parameter Real eta = 1 "Value of efficiency";
        parameter SI.Length P(displayUnit = "mm") "Pitch of spindle";
        parameter SI.Force C_am = 50000 "Dynamic load rating";
        parameter SI.Force F_pr "Preload force";
        Modelica.Mechanics.Rotational.Components.IdealGearR2T idealGearR2T(ratio = 2 * Modelica.Constants.pi / P, useSupportR = false, useSupportT = true) annotation(Placement(transformation(extent = {{20, -10}, {40, 10}})));
        RotationalComponents.BasicComponents.EfficiencyFactor efficiencyFactor(eta = eta) annotation(Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{60, -10}, {80, 10}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {0, -30})));
        NominalLifeRating nominalLifeRating(C_am = C_am, F_pr = F_pr, temp_F(fixed = true), temp_w(fixed = true)) annotation(Placement(transformation(extent = {{80, -60}, {100, -40}})));
      public
        Modelica.Blocks.Math.UnitConversions.To_rpm speedSensor_rpm annotation(Placement(transformation(extent = {{20, -90}, {40, -70}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flangeT
          "(right) driven flange (flange axis directed out of cut plane)"                                                            annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flangeR
          "Left flange of shaft"                                                         annotation(Placement(transformation(extent = {{-114, -10}, {-94, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Support supportT1
          "Translational support/housing of component"                                                             annotation(Placement(transformation(extent = {{-112, -60}, {-92, -40}}), iconTransformation(extent = {{-112, -60}, {-92, -40}})));
        HelpBlocks.AbsMax speed_absMax
          annotation (Placement(transformation(extent={{56,-90},{76,-70}})));
        HelpBlocks.AbsMax force_absMax
          annotation (Placement(transformation(extent={{80,-32},{100,-12}})));
      equation
        connect(idealGearR2T.flangeT, forceSensor.flange_a) annotation(Line(points = {{40, 0}, {60, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, idealGearR2T.flangeR) annotation(Line(points = {{1.77636e-015, -20}, {1.77636e-015, 0}, {20, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(forceSensor.f, nominalLifeRating.u1) annotation(Line(points = {{62, -11}, {62, -44}, {78, -44}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, nominalLifeRating.u2) annotation(Line(points = {{-1.9984e-015, -41}, {-1.9984e-015, -56}, {78, -56}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, speedSensor_rpm.u) annotation(Line(points = {{-1.9984e-015, -41}, {-1.9984e-015, -80}, {18, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(efficiencyFactor.flange_b, idealGearR2T.flangeR) annotation(Line(points = {{-50, 0}, {20, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(forceSensor.flange_b, flangeT) annotation(Line(points = {{80, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(efficiencyFactor.flange_a, flangeR) annotation(Line(points = {{-70, 0}, {-104, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGearR2T.supportT, supportT1) annotation(Line(points = {{40, -10}, {40, -50}, {-102, -50}}, color = {0, 0, 0}, pattern = LinePattern.None, smooth = Smooth.None));
        connect(speedSensor_rpm.y, speed_absMax.u) annotation (Line(
            points={{41,-80},{54,-80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force_absMax.u, nominalLifeRating.u1) annotation (Line(
            points={{78,-22},{62,-22},{62,-44},{78,-44}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-60, 26}, {60, 20}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Rectangle(extent=  {{-100, 20}, {80, -20}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-60, 0}, {60, 0}, {60, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-60, 20}, {-40, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-40, 20}, {-20, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-20, 20}, {0, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{0, 20}, {20, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{20, 20}, {40, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{40, 20}, {60, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Rectangle(extent=  {{-60, 0}, {60, -26}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Rectangle(extent=  {{-40, 40}, {-20, -40}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-60, 0}, {-100, 0}, {-102, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-60, 20}, {-60, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{60, 20}, {60, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-80, 20}, {-60, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-100, 20}, {-80, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{60, 20}, {80, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{60, 26}, {60, 26}, {100, 26}, {100, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-100, -50}, {-74, -50}, {-74, -20}}, color=  {0, 0, 0}, smooth=  Smooth.None), Text(extent=  {{-150, 100}, {150, 60}}, lineColor=  {0, 0, 255}, textString=  "%name")}));
      end BallScrewNut;
    end BallScrew;

    package RackPinion
      model RackPinion "Model of a rack-pinion drive"
        parameter Integer mod "Module" annotation(Dialog(group = "Engineering data"));
        parameter Integer z "Number of teeth" annotation(Dialog(group = "Engineering data"));
        parameter SI.Diameter d "Diameter of gearweheel" annotation(Dialog(group = "Geometric data"));
        parameter SI.Diameter d_1 "Bore diameter of gearwheel" annotation(Dialog(group = "Geometric data"));
        parameter SI.Mass m "Mass of gearwheel" annotation(Dialog(group = "Physical data"));
        parameter SI.MomentOfInertia J = m * (d ^ 2 + d_1 ^ 2) / 8
          "Moment of inertia of gearwheel"                                                          annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter Real eta = 0.9 "Value of efficiency" annotation(Dialog(group = "Physical data"));
        parameter Real K_A = 1.5 "Load factor (see e.g. Atlanta)" annotation(Dialog(group = "Engineering data"));
        parameter Real S_B = 1.2 "Safety coefficient (see e.g. Atlanta)" annotation(Dialog(group = "Engineering data"));
        parameter Real f_n = 1.0 "Life-time factor (see e.g. Atlanta)" annotation(Dialog(group = "Engineering data"));
        parameter Real L_KHb = 1.1
          "Linear oad distribution factor (see e.g. Atlanta)"                          annotation(Dialog(group = "Engineering data"));
        parameter SI.Force F_perm_0 = 1000
          "Maximum feed force (Catalogue value)"                                  annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_perm = F_perm_0 / (K_A * S_B * f_n * L_KHb)
          "Permissible force (with factors)"                                                                annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        SI.Force con_F_ax "Constraint for axial force";
        TranslationalComponents.EfficiencyFactor_Translational eff(eta = eta) annotation(Placement(transformation(extent = {{20, -10}, {40, 10}})));
        Modelica.Mechanics.Rotational.Components.IdealGearR2T idealGearR2T(ratio = 2 / d) annotation(Placement(transformation(extent = {{-40, -10}, {-20, 10}})));
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J = J) annotation(Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{-8, -10}, {12, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_r
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-90, -10}, {-70, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax annotation(Placement(transformation(extent = {{50, -10}, {70, 10}})));
        HelpBlocks.AbsMax absMax
          annotation (Placement(transformation(extent={{2,-46},{22,-26}})));
      equation
        con_F_ax = absMax.y - F_perm;
        connect(inertia.flange_b, idealGearR2T.flangeR) annotation(Line(points = {{-50, 0}, {-40, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGearR2T.flangeT, forceSensor.flange_a) annotation(Line(points = {{-20, 0}, {-8, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.flange_b, eff.flange_ax_a) annotation(Line(points = {{12, 0}, {21, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(inertia.flange_a, flange_r) annotation(Line(points = {{-70, 0}, {-80, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(eff.flange_ax_b, flange_ax) annotation(Line(points = {{41, 0}, {60, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.f, absMax.u) annotation (Line(
            points={{-6,-11},{-6,-36},{0,-36}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-44, 6}, {-2, -36}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-50, 12}, {4, -42}}, lineColor=  {0, 0, 0}), Rectangle(extent=  {{-80, 20}, {80, 6}}, lineColor=  {0, 0, 0}), Line(points=  {{-80, 12}, {60, 12}, {80, 12}}, color=  {0, 0, 0}, smooth=  Smooth.None, pattern=  LinePattern.Dash), Line(points=  {{-80, -14}, {-80, 0}, {-98, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{90, 0}, {84, 0}, {84, 16}, {80, 16}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-24, -14}, {-78, -14}, {-80, -14}}, color=  {0, 0, 0}, smooth=  Smooth.None)}));
      end RackPinion;
    end RackPinion;

    package ToothedBeltAxis
      model ToothedBeltAxis
        "Model of a toothed belt axis from suppliers such as Festo"
        parameter SI.Diameter d(displayUnit = "mm") = 0.01846
          "Diameter of wheel"                                                     annotation(Dialog(group = "Geometric data"));
        parameter SI.Diameter l = 1 "Working stroke" annotation(Dialog(group = "Geometric data"));
        parameter Integer n_slides = 1 "Number of slides" annotation(Dialog(group = "Physical data"));
        parameter SI.MomentOfInertia J_0 = 16.94e-6
          "Moment of inertia of wheels"                                           annotation(Dialog(group = "Physical data (Inertia)"));
        parameter Real J_H(unit = "kg.m") = 2.6e-6
          "Moment of inertia per metre stroke"                                          annotation(Dialog(group = "Physical data (Inertia)"));
        parameter SI.MomentOfInertia J_W = 3.56e-6
          "Moment of inertia of each additional slides"                                          annotation(Dialog(group = "Physical data (Inertia)"));
        parameter SI.MomentOfInertia J_F = 0
          "Moment of inertia of clamping unit"                                    annotation(Dialog(group = "Physical data (Inertia)"));
        parameter SI.MomentOfInertia J = J_0 + n_slides * J_W + l * J_H + J_F
          "Moment of inertia of axis"                                                                     annotation(Dialog(group = "Physical data", tab = "Calculated"));
        parameter SI.Force F_Fr = 8 "No load resistance to shifting" annotation(Dialog(group = "Physical data"));
        parameter SI.Force F_perm = 50 "Maximum feed force (Catalogue value)" annotation(Dialog(group = "Limiting data"));
        parameter SI.Acceleration a_perm = 50
          "Maximum acceleration (Catalogue value)"                                     annotation(Dialog(group = "Limiting data"));
        parameter SI.Velocity v_perm = 3 "Maximum velocity (Catalogue value)" annotation(Dialog(group = "Limiting data"));
        parameter Real expan = 9.4e-4
          "Expansion at max. feed force of toothed belt"                             annotation(Dialog(group = "Physical data"));
        parameter SI.TranslationalSpringConstant c_ax = F_perm / (expan * l)
          "Rigidity of toothed belt"                                                                    annotation(Dialog(group = "Physical data"));
        parameter SI.Mass m_ref = 1000 "Reference mass for eigenfrequency" annotation(Dialog(group = "Engineering data"));
        parameter SI.Frequency f_bs = 1 / (2 * Constants.pi) * sqrt(c_ax / m_ref)
          "Eigenfrequncy of axis"                                                                         annotation(Dialog(group = "Engineering data", tab = "Calculated"));
        parameter SI.Frequency f_min = 10
          "Minimum eigenfrequncy of axis. Constraint f_bs > f_min"                                 annotation(Dialog(group = "Limiting data"));
        SI.Force con_F_ax "Constraint for axial force";
        SI.Acceleration con_a "Constraint for acceleration";
        SI.Velocity con_v "Constraint for velocity";
        SI.Frequency con_f "Constraint for eigenfrequncy of axis";
        Modelica.Mechanics.Rotational.Components.IdealGearR2T idealGearR2T(ratio = 2 / d) annotation(Placement(transformation(extent = {{-48, -10}, {-28, 10}})));
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J = J) annotation(Placement(transformation(extent = {{-78, -10}, {-58, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_r
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
        TranslationalComponents.Sensors.TranslationalSensor s annotation(Placement(transformation(extent = {{-6, -10}, {14, 10}})));
        TranslationalComponents.NoLoad_Friction Friction(F_Fr_const = F_Fr) annotation(Placement(transformation(extent = {{38, -10}, {58, 10}})));
      equation
        con_F_ax = s.f_max - F_perm;
        con_a = s.a_max - a_perm;
        con_v = s.v_max - v_perm;
        con_f = f_min - f_bs;
        connect(inertia.flange_b, idealGearR2T.flangeR) annotation(Line(points = {{-58, 0}, {-48, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(inertia.flange_a, flange_r) annotation(Line(points = {{-78, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGearR2T.flangeT, s.flange_a) annotation(Line(points = {{-28, 0}, {-6.2, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(s.flange_b, Friction.flange_ax_a) annotation(Line(points = {{14, 0}, {38, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(Friction.flange_ax_b, flange_ax) annotation(Line(points = {{58, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-74, 24}, {-26, -24}}, lineColor=  {0, 0, 0}), Ellipse(extent=  {{28, 24}, {76, -24}}, lineColor=  {0, 0, 0}), Line(points=  {{-50, 24}, {52, 24}}, color=  {0, 0, 0}, smooth=  Smooth.None, pattern=  LinePattern.Dot), Line(points=  {{-50, -24}, {52, -24}}, color=  {0, 0, 0}, smooth=  Smooth.None, pattern=  LinePattern.Dot), Line(points=  {{-50, 2}, {-50, -2}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-48, 0}, {-52, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{52, 2}, {52, -2}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{54, 0}, {50, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Rectangle(extent=  {{-16, 30}, {12, 24}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 28}, {100, 28}, {100, 6}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-50, 0}, {-98, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None)}));
      end ToothedBeltAxis;
    end ToothedBeltAxis;
  end LinearActuators;

  package TranslationalComponents


    model Force_PT1
      "External force acting on a drive train element as input signal"
      extends
        Modelica.Mechanics.Translational.Interfaces.PartialElementaryOneFlangeAndSupport2;
      Modelica.Blocks.Interfaces.RealInput f "Driving force as input signal" annotation(Placement(transformation(extent = {{-140, -20}, {-100, 20}}, rotation = 0)));
      parameter SI.Time T = 0.001 "Time for PT1";
    public
      Modelica.Blocks.Continuous.FirstOrder PT1(k = 1, T = T) annotation(Placement(transformation(extent = {{-80, -10}, {-60, 10}})));
    equation
      flange.f = -PT1.y;
      connect(f, PT1.u) annotation(Line(points = {{-120, 0}, {-82, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation(Documentation(info = "<html>
<p>
The input signal \"f\" in [N] characterizes an <i>external
force</i> which acts (with positive sign) at a flange,
i.e., the component connected to the flange is driven by force f.
</p>
<p>
Input signal f can be provided from one of the signal generator
blocks of Modelica.Blocks.Source.
</p>

</HTML>
"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Polygon(points = {{-100, 10}, {20, 10}, {20, 41}, {90, 0}, {20, -41}, {20, -10}, {-100, -10}, {-100, 10}}, lineColor = {0, 127, 0}, fillColor = {215, 215, 215},
                fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{-150, -32}, {-80, -62}}, lineColor = {0, 0, 0}, textString = "f"), Text(extent = {{-150, 90}, {150, 50}}, textString = "%name", lineColor = {0, 0, 255}), Line(points = {{-30, -60}, {30, -60}}, color = {0, 0, 0}), Line(points = {{0, -60}, {0, -101}}, color = {0, 0, 0}), Line(points = {{-30, -80}, {-10, -60}}, color = {0, 0, 0}), Line(points = {{-10, -80}, {10, -60}}, color = {0, 0, 0}), Line(points = {{10, -80}, {30, -60}}, color = {0, 0, 0}), Polygon(points = {{-61, -50}, {-30, -40}, {-30, -60}, {-61, -50}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-31, -50}, {50, -50}}, color = {0, 0, 0}), Line(points = {{-50, -80}, {-30, -60}}, color = {0, 0, 0})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics));
    end Force_PT1;

    model LinearGuide
      "Model of a linear guide friction and life-time calculation"
      parameter SI.Mass m(min = 0, start = 1) = 1 "Mass of the sliders" annotation(Dialog(group = "Physical data"));
      parameter SI.CoefficientOfFriction mu = 0.004 "Friction coefficient" annotation(Dialog(group = "Physical data"));
      parameter SI.Force F_guideCover = 500 "Friction force of guide cover" annotation(Dialog(group = "Physical data"));
      parameter SI.Velocity v_ref = 0.001
        "Velocity after which the friction force is approx. equal to mu*F_N"                                   annotation(Dialog(tab = "Advanced"));
      parameter SI.Force C_100(displayUnit = "kN") = 28000
        "Dynamic load rating for a nominal lifetime 100 km"                                                    annotation(Dialog(group = "Limiting data"));
      parameter SI.Force C_0(displayUnit = "kN") = 65000 "Static load rating" annotation(Dialog(group = "Limiting data"));
      parameter Real S_0 = 1 "Safety factor for static load rating" annotation(Dialog(group = "Limiting data"));
      parameter SI.Velocity v_lim = 4 "Maximum velocity" annotation(Dialog(group = "Limiting data"));
      parameter SI.Acceleration a_lim = 100 "Maximum acceleration" annotation(Dialog(group = "Limiting data"));
      parameter Sensors.RailSystemType railSystemType "Type of system" annotation(Dialog(group = "Engineering data"));
      parameter SI.Conversions.NonSIunits.Time_hour L_h_min = 10000
        "Minimum life expectancy of nut"                                                             annotation(Dialog(group = "Limiting data"));
      Interfaces.Flange_b flange_b annotation(Placement(transformation(extent = {{70, -10}, {90, 10}})));
      Interfaces.FlangeChangeRight flangeChangeRight annotation(Placement(transformation(extent = {{16, 20}, {36, 0}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(Placement(transformation(extent = {{-90, -10}, {-70, 10}}), iconTransformation(extent = {{-90, -10}, {-70, 10}})));
      SimpleFriction rollerFriction(mu = mu, v_ref = v_ref) annotation(Placement(transformation(extent = {{-34, -10}, {-14, 10}})));
      NoLoad_Friction coverFriction(v_ref = v_ref, F_Fr_const = F_guideCover) annotation(Placement(transformation(extent = {{-64, -10}, {-44, 10}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_N annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-24, -20}), iconTransformation(extent = {{-10, -66}, {10, -46}})));
      Sensors.LinearRollerGuideSensor linGuideReq(p = p, C_100 = C_100) annotation(Placement(transformation(extent = {{46, -10}, {66, 10}})));
      Modelica.Mechanics.Translational.Components.Mass guideMass(m = m, s(start = 0), v(start = 0), a(start = 0)) annotation(Placement(transformation(extent = {{-6, -10}, {14, 10}})));
      SI.Velocity con_v;
      SI.Acceleration con_a;
      SI.Conversions.NonSIunits.Time_hour con_L_h;
      SI.Force con_C_0;
    protected
      parameter Real p = if railSystemType == Sensors.RailSystemType.RollerRailSystem then 10 / 3 else 3
        "Life expectancy exponent";
    equation
      con_v = linGuideReq.v_max - v_lim;
      con_a = linGuideReq.a_max - a_lim;
      con_C_0 = S_0 * linGuideReq.P_max - C_0;
      con_L_h = linGuideReq.L_h - L_h_min;
      connect(flangeChangeRight.flange_N, rollerFriction.flange_N_a) annotation(Line(points = {{21, 20}, {-24, 20}, {-24, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(coverFriction.flange_ax_b, rollerFriction.flange_ax_a) annotation(Line(points = {{-44, 0}, {-34, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(coverFriction.flange_ax_a, flange_a) annotation(Line(points = {{-64, 0}, {-80, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(rollerFriction.flange_N_b, flange_N) annotation(Line(points = {{-24, -10}, {-24, -20}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(linGuideReq.flange_b1, flange_b) annotation(Line(points = {{66, 0}, {80, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(rollerFriction.flange_ax_b, guideMass.flange_a) annotation(Line(points = {{-14, 0}, {-6, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(guideMass.flange_b, flangeChangeRight.flange_ax) annotation(Line(points = {{14, 0}, {21, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChangeRight.flange_b, linGuideReq.flange_a1) annotation(Line(points = {{31, 10}, {38, 10}, {38, 0}, {46, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, -6}, {80, -46}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-60, 14}, {6, -22}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward), Line(points = {{-80, -22}, {80, -22}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{16, 16}, {68, 16}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{-46, 54}, {-46, 20}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{-90, 0}, {-60, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{6, 0}, {92, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-150, 100}, {150, 60}}, textString = "%name", lineColor = {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
    end LinearGuide;

    model Slide "Calculates normal and tangential force"
      parameter SI.Angle alpha = 0 "Angle of Inclination";
      parameter SI.Mass m = 1000 "Moved mass";
      parameter SI.Force F_T = m * Constants.g_n * sin(alpha);
      parameter SI.Force F_N = m * Constants.g_n * cos(alpha);
      Modelica.Mechanics.Translational.Components.Mass mass(m = m, s(start = 0), v(start = 0), a(start = 0)) annotation(Placement(transformation(extent = {{-2, 0}, {18, 20}})));
      Modelica.Mechanics.Translational.Sources.ConstantForce constantForce(f_constant = F_N) annotation(Placement(transformation(extent = {{-20, -50}, {0, -30}})));
      Modelica.Mechanics.Translational.Sources.ConstantForce force_T(f_constant = -F_T) annotation(Placement(transformation(extent = {{-80, 14}, {-60, 34}})));
      Interfaces.Flange_a flange_a annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
      Interfaces.FlangeChangeLeft flangeChange annotation(Placement(transformation(extent = {{-66, -10}, {-46, 10}})));
      Interfaces.Flange_b flange_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
      Interfaces.FlangeChangeRight flangeChangeRight annotation(Placement(transformation(extent = {{62, -10}, {82, 10}})));
    equation
      connect(flangeChange.flange_a, flange_a) annotation(Line(points = {{-61, 0}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChange.flange_ax, mass.flange_a) annotation(Line(points = {{-51, 10}, {-2, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(mass.flange_b, flangeChangeRight.flange_ax) annotation(Line(points = {{18, 10}, {67, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(force_T.flange, mass.flange_a) annotation(Line(points = {{-60, 24}, {-2, 24}, {-2, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChange.flange_N, flangeChangeRight.flange_N) annotation(Line(points = {{-51, -10}, {7.5, -10}, {7.5, -10}, {67, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChangeRight.flange_b, flange_b) annotation(Line(points = {{77, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(constantForce.flange, flangeChangeRight.flange_N) annotation(Line(points = {{0, -40}, {18, -40}, {18, -10}, {67, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-60, 40}, {60, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{-60, 40}, {60, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid, textString = "m"), Line(points = {{-60, 0}, {-98, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{90, 0}, {60, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-150, 100}, {150, 60}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-160, -50}, {160, -86}}, lineColor = {0, 0, 0}, textString = "m=%m")}));
    end Slide;

    model SimpleFriction
      "Simple friction model with atan-function. No static friction force is assumed"
      parameter SI.CoefficientOfFriction mu = 0.004 "Friction coefficient";
      SI.Velocity v "Velocity of guide";
      SI.Force F_Fr "Friction force";
      SI.Force F_N "Normal force";
      parameter SI.Velocity v_ref = 0.001
        "Velocity after which the friction force is approx. equal to mu*F_N";
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_ax_a annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_N_a annotation(Placement(transformation(extent = {{-10, 90}, {10, 110}}), iconTransformation(extent = {{-10, 90}, {10, 110}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_N_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{-10, -110}, {10, -90}})));
    equation
      flange_ax_a.s = flange_ax_b.s;
      flange_N_a.s = flange_N_b.s;
      flange_N_a.f + flange_N_b.f = 0;
      v = der(flange_ax_a.s);
      F_N = max(0, flange_N_a.f);
      F_Fr = -2 / Constants.pi * F_N * mu * atan(v / v_ref);
      flange_ax_a.f + flange_ax_b.f + F_Fr = 0;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, -22}, {80, -46}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{-32, 14}, {34, -22}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward), Line(points = {{-80, -22}, {80, -22}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{34, 4}, {86, 4}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{40, 38}, {80, -2}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward, textString = "v"), Line(points = {{-22, 54}, {-22, 20}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{-94, 56}, {-24, 26}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward, textString = "F_N"), Line(points = {{-90, 0}, {-32, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{6, 0}, {92, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-148, -44}, {152, -84}}, lineColor = {0, 0, 255}, textString = "%name"), Line(points = {{0, 90}, {0, 12}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{0, -90}, {0, -22}}, color = {0, 0, 0}, smooth = Smooth.None)}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
    end SimpleFriction;

    model NoLoad_Friction
      SI.Velocity v "Velocity of guide";
      parameter SI.Force F_Fr_const "Friction force";
      SI.Force F_Fr "Directed friction force";
      parameter SI.Velocity v_ref = 0.001
        "Velocity after which the friction force is approx. equal to mu*F_N";
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_ax_a annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
    equation
      flange_ax_a.s = flange_ax_b.s;
      v = der(flange_ax_a.s);
      F_Fr = -2 / Constants.pi * F_Fr_const * atan(v / v_ref);
      flange_ax_a.f + flange_ax_b.f + F_Fr = 0;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, -22}, {80, -46}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{-32, 14}, {34, -22}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward), Line(points = {{-80, -22}, {80, -22}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{34, 4}, {86, 4}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{40, 38}, {80, -2}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward, textString = "v"), Line(points = {{16, 20}, {-22, 20}}, color = {0, 0, 0}, smooth = Smooth.None, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{-36, 60}, {34, 30}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward, textString = "F_F"), Line(points = {{-90, 0}, {-32, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{6, 0}, {92, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Text(extent = {{-148, -44}, {152, -84}}, textString = "%name", lineColor = {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
    end NoLoad_Friction;



    model EfficiencyFactor_Translational
      "Calculates the friction roce by using the efficiency factor. Atan function is used to obtain a continous model"
      parameter Real eta = 1 "efficiency factor";
      SI.Velocity v "Velocity of flange";
      SI.Force F_F(start = 0) "Friction force";
      parameter SI.AngularVelocity v_ref = 0.001
        "Velocity after which the friction torque can is approx. equal to (1-eta)*flange_a.tau";
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_ax_a annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-100, -10}, {-80, 10}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
    equation
      flange_ax_a.s = flange_ax_b.s;
      v = der(flange_ax_a.s);
      F_F = -2 / Constants.pi * (1 - eta) * sign(flange_ax_a.f) * flange_ax_a.f * atan(v / v_ref);
      flange_ax_a.f + flange_ax_b.f + F_F = 0;
      annotation(Icon(graphics={  Text(extent = {{-150, 100}, {150, 60}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-144, -50}, {156, -90}}, lineColor = {0, 0, 0}, textString = "eta=%eta"), Polygon(points = {{-99, 50}, {-70, 50}, {-70, 90}, {-80, 90}, {-60, 110}, {-40, 90}, {-50, 90}, {-50, 30}, {-99, 30}, {-99, 50}}, lineColor = {0, 0, 0}, fillColor = {255, 0, 0},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-70, -22}, {90, -46}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Backward), Rectangle(extent = {{-22, 14}, {44, -22}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Backward), Line(points = {{-70, -22}, {90, -22}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-80, 0}, {-22, 0}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{16, 0}, {102, 0}}, color = {0, 0, 0}, smooth = Smooth.None)}));
    end EfficiencyFactor_Translational;

    model SpringDamper "Spring damper behavior for horizontal force"
      parameter SI.TranslationalSpringConstant c(final min = 0, start = 1)
        "Spring constant";
      parameter SI.TranslationalDampingConstant d(final min = 0, start = 1)
        "Damping constant";
      parameter SI.Position s_rel0 = 0 "Unstretched spring length";
      Modelica.Mechanics.Translational.Components.SpringDamper springDamper(c = c, d = d, s_rel0 = s_rel0) annotation(Placement(transformation(extent = {{-10, 0}, {10, 20}})));
      Interfaces.Flange_a flange_a annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Interfaces.Flange_b flange_b annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
      Interfaces.FlangeChangeLeft flangeChangeLeft annotation(Placement(transformation(extent = {{-66, -10}, {-46, 10}})));
      Interfaces.FlangeChangeRight flangeChangeRight annotation(Placement(transformation(extent = {{54, -10}, {74, 10}})));
    equation
      connect(flangeChangeLeft.flange_a, flange_a) annotation(Line(points = {{-61, 0}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChangeLeft.flange_ax, springDamper.flange_a) annotation(Line(points = {{-51, 10}, {-10, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(springDamper.flange_b, flangeChangeRight.flange_ax) annotation(Line(points = {{10, 10}, {59, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChangeLeft.flange_N, flangeChangeRight.flange_N) annotation(Line(points = {{-51, -10}, {59, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(flangeChangeRight.flange_b, flange_b) annotation(Line(points = {{69, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      annotation(Documentation(info = "<html>
<p>
A <i>spring and damper element connected in parallel</i>.
The component can be
connected either between two sliding masses to describe the elasticity
and damping, or between a sliding mass and the housing (model Fixed),
to describe a coupling of the sliding mass with the housing via a spring/damper.
</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points = {{-80, 40}, {-60, 40}, {-45, 10}, {-15, 70}, {15, 10}, {45, 70}, {60, 40}, {80, 40}}, color = {0, 0, 0}), Line(points = {{-80, 40}, {-80, -70}}, color = {0, 0, 0}), Line(points = {{-80, -70}, {-52, -70}}, color = {0, 0, 0}), Rectangle(extent = {{-52, -49}, {38, -91}}, lineColor = {0, 0, 0}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-52, -49}, {68, -49}}, color = {0, 0, 0}), Line(points = {{-51, -91}, {69, -91}}, color = {0, 0, 0}), Line(points = {{38, -70}, {80, -70}}, color = {0, 0, 0}), Line(points = {{80, 40}, {80, -70}}, color = {0, 0, 0}), Line(points = {{-90, 0}, {-80, 0}}, color = {0, 0, 0}), Line(points = {{80, 0}, {90, 0}}, color = {0, 0, 0}), Polygon(points = {{53, -18}, {23, -8}, {23, -28}, {53, -18}}, lineColor = {128, 128, 128}, fillColor = {128, 128, 128},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-57, -18}, {23, -18}}, color = {0, 0, 0}), Text(extent = {{-150, 120}, {150, 80}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-150, -135}, {150, -165}}, lineColor = {0, 0, 0}, textString = "d=%d"), Text(extent = {{-150, -100}, {150, -130}}, lineColor = {0, 0, 0}, textString = "c=%c"), Line(visible = useHeatPort, points = {{-100, -100}, {-100, -80}, {-5, -80}}, color = {191, 0, 0}, pattern = LinePattern.Dot, smooth = Smooth.None)}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
    end SpringDamper;

    package Interfaces
      connector Flange_a
        "(left) 2D translational flange (flange axis directed INTO cut plane, e. g. from left to right)"
        SI.Position s_ax "Axial position of flange";
        flow SI.Force f_ax "Axial cut force directed into flange";
        SI.Position s_N "Normal position of flange";
        flow SI.Force f_N "Normal cut force directed into flange";
        annotation(defaultComponentName = "flange_a", Documentation(info = "<html>
This is a flange for 1D translational mechanical systems. In the cut plane of
the flange a unit vector n, called flange axis, is defined which is directed
INTO the cut plane, i. e. from left to right. All vectors in the cut plane are
resolved with respect to
this unit vector. E.g. force f characterizes a vector which is directed in
the direction of n with value equal to f. When this flange is connected to
other 1D translational flanges, this means that the axes vectors of the connected
flanges are identical.
</p>
<p>
The following variables are transported through this connector:
<pre>
  s: Absolute position of the flange in [m]. A positive translation
     means that the flange is translated along the flange axis.
  f: Cut-force in direction of the flange axis in [N].
</pre>
</HTML>
"),       Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,127,0},
                      fillColor={0,127,0},
                      fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={Rectangle(
                      extent={{-40,-40},{40,40}},
                      lineColor={0,127,0},
                      fillColor={0,127,0},
                      fillPattern=FillPattern.Solid),Text(
                      extent={{-160,110},{40,50}},
                      lineColor={0,127,0},
                      textString="%name")}));
      end Flange_a;

      connector Flange_b
        "(right) 2D translational flange (flange axis directed OUT OF cut plane)"
        SI.Position s_ax "Axial position of flange";
        flow SI.Force f_ax "Axial cut force directed into flange";
        SI.Position s_N "Normal position of flange";
        flow SI.Force f_N "Normal cut force directed into flange";
        annotation(defaultComponentName = "flange_b", Documentation(info = "<html>
This is a flange for 1D translational mechanical systems. In the cut plane of
the flange a unit vector n, called flange axis, is defined which is directed
OUT OF the cut plane. All vectors in the cut plane are resolved with respect to
this unit vector. E.g. force f characterizes a vector which is directed in
the direction of n with value equal to f. When this flange is connected to
other 1D translational flanges, this means that the axes vectors of the connected
flanges are identical.
</p>
<p>
The following variables are transported through this connector:
<pre>
  s: Absolute position of the flange in [m]. A positive translation
     means that the flange is translated along the flange axis.
  f: Cut-force in direction of the flange axis in [N].
</pre>
</HTML>
"),       Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Rectangle(
                extent={{-100,-100},{100,100}},
                lineColor={0,127,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Rectangle(
                extent={{-40,-40},{40,40}},
                lineColor={0,127,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-40,110},{160,50}},
                lineColor={0,127,0},
                textString="%name")}));
      end Flange_b;

      model FlangeChangeLeft
        Flange_a flange_a annotation(Placement(transformation(extent = {{-60, -10}, {-40, 10}}), iconTransformation(extent = {{-60, -10}, {-40, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_ax annotation(Placement(transformation(extent = {{40, 90}, {60, 110}}), iconTransformation(extent = {{40, 90}, {60, 110}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_N annotation(Placement(transformation(extent = {{40, -110}, {60, -90}}), iconTransformation(extent = {{40, -110}, {60, -90}})));
      equation
        flange_a.s_ax = flange_ax.s;
        flange_a.f_ax + flange_ax.f = 0;
        flange_a.s_N = flange_N.s;
        flange_a.f_N + flange_N.f = 0;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-40, 110}, {40, -110}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{62, 116}, {-62, 68}}, lineColor=  {0, 0, 0},
                  lineThickness=                                                                                                    1, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "ax"), Text(extent=  {{66, -66}, {-60, -112}}, lineColor=  {0, 0, 0},
                  lineThickness=                                                                                                    1, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "N")}));
      end FlangeChangeLeft;

      model FlangeChangeRight
        Flange_b flange_b annotation(Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_ax annotation(Placement(transformation(extent = {{-60, 90}, {-40, 110}}), iconTransformation(extent = {{-60, 90}, {-40, 110}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_N annotation(Placement(transformation(extent = {{-60, 92}, {-40, 112}}), iconTransformation(extent = {{-60, -110}, {-40, -90}})));
      equation
        flange_b.s_ax = flange_ax.s;
        flange_b.f_ax + flange_ax.f = 0;
        flange_b.s_N = flange_N.s;
        flange_b.f_N + flange_N.f = 0;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-40, 110}, {40, -110}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{62, 116}, {-62, 68}}, lineColor=  {0, 0, 0},
                  lineThickness=                                                                                                    1, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "Ax"), Text(extent=  {{66, -64}, {-58, -112}}, lineColor=  {0, 0, 0},
                  lineThickness=                                                                                                    1, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, textString=  "N")}));
      end FlangeChangeRight;
    end Interfaces;

    package Sensors
      model TranslationalSensor
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}})));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-60, 30}, {-40, 50}})));
        Modelica.Mechanics.Translational.Sensors.AccSensor accSensor annotation(Placement(transformation(extent = {{-60, 70}, {-40, 90}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a
          "(left) driving flange (flange axis directed in to cut plane, e. g. from left to right)"
                                                                                                              annotation(Placement(transformation(extent = {{-112, -10}, {-92, 10}}), iconTransformation(extent = {{-112, -10}, {-92, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b
          "(right) driven flange (flange axis directed out of cut plane)"                                                             annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Interfaces.RealOutput a_actual
          "Absolute acceleration of flange as output signal"                                              annotation(Placement(transformation(extent = {{100, 50}, {120, 70}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {50, -110})));
        Modelica.Blocks.Interfaces.RealOutput a_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, 70}, {120, 90}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {80, -110})));
        Modelica.Blocks.Interfaces.RealOutput v_actual
          "Absolute velocity of flange as output signal"                                              annotation(Placement(transformation(extent = {{100, 10}, {120, 30}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-20, -110})));
        Modelica.Blocks.Interfaces.RealOutput v_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{46, -114}, {66, -94}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {20, -110})));
        Modelica.Blocks.Interfaces.RealOutput f_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, -40}, {120, -20}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-50, -110})));
        Modelica.Blocks.Interfaces.RealOutput f_actual
          "Connector of Real output signal"                                              annotation(Placement(transformation(extent = {{100, -70}, {120, -50}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        HelpBlocks.AbsMax absMax3
          annotation (Placement(transformation(extent={{0,70},{20,90}})));
        HelpBlocks.AbsMax absMax
          annotation (Placement(transformation(extent={{0,30},{20,50}})));
        HelpBlocks.AbsMax absMax1
          annotation (Placement(transformation(extent={{10,-40},{30,-20}})));
      equation
        connect(forceSensor.flange_a, flange_a) annotation(Line(points = {{-10, 0}, {-102, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.flange_b, flange_b) annotation(Line(points = {{10, 0}, {54, 0}, {54, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flange_a) annotation(Line(points = {{-60, 40}, {-60, 0}, {-102, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(accSensor.flange, flange_a) annotation(Line(points = {{-60, 80}, {-60, 0}, {-102, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(accSensor.a, a_actual) annotation(Line(points = {{-39, 80}, {-36, 80}, {-36, 60}, {110, 60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.v, v_actual) annotation(Line(points = {{-39, 40}, {-36, 40}, {-36, 20}, {110, 20}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(forceSensor.f, f_actual) annotation(Line(points = {{-8, -11}, {-8, -60}, {110, -60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMax3.y, a_max) annotation (Line(
            points={{21,80},{110,80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax3.u, a_actual) annotation (Line(
            points={{-2,80},{-36,80},{-36,60},{110,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax.u, speedSensor.v) annotation (Line(
            points={{-2,40},{-39,40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax.y, v_max) annotation (Line(
            points={{21,40},{36,40},{36,-104},{56,-104}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax1.y, f_max) annotation (Line(
            points={{31,-30},{110,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax1.u, f_actual) annotation (Line(
            points={{8,-30},{-8,-30},{-8,-60},{110,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-74, -60}, {66, 20}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Polygon(points=  {{-4, -40}, {-14, -16}, {6, -16}, {-4, -40}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-4, 0}, {-4, -16}}, color=  {0, 0, 0}), Line(points=  {{-74, 0}, {-4, 0}}, color=  {0, 0, 0}), Line(points=  {{-54, -40}, {-54, -60}}, color=  {0, 0, 0}), Line(points=  {{-34, -40}, {-34, -60}}, color=  {0, 0, 0}), Line(points=  {{-14, -40}, {-14, -60}}, color=  {0, 0, 0}), Line(points=  {{6, -40}, {6, -60}}, color=  {0, 0, 0}), Line(points=  {{26, -40}, {26, -60}}, color=  {0, 0, 0}), Line(points=  {{46, -40}, {46, -60}}, color=  {0, 0, 0}), Line(points=  {{-90, -80}, {-10, -80}}, color=  {0, 0, 0}), Polygon(points=  {{16, -82}, {-14, -72}, {-14, -92}, {16, -82}}, lineColor=  {128, 128, 128}, fillColor=  {128, 128, 128},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-74, 0}, {-94, 0}}, color=  {0, 0, 0}), Line(points=  {{66, 0}, {96, 0}}, color=  {0, 0, 127}), Text(extent=  {{-154, 100}, {146, 60}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{-74, 0}, {-94, 0}}, color=  {0, 0, 0})}));
      end TranslationalSensor;

      model LinearRollerGuideSensor
        parameter Real p "Life expectancy exponent";
        parameter SI.Force C_100
          "Dynmamic load rating for a nominal lifetime 100 km";
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{-8, -20}, {12, 0}})));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-60, 30}, {-40, 50}})));
        Modelica.Mechanics.Translational.Sensors.AccSensor accSensor annotation(Placement(transformation(extent = {{-60, 70}, {-40, 90}})));
        Modelica.Blocks.Interfaces.RealOutput a_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, 70}, {120, 90}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {20, -110})));
        Modelica.Blocks.Interfaces.RealOutput v_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, 30}, {120, 50}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-20, -110})));
        Interfaces.FlangeChangeLeft flangeChangeLeft annotation(Placement(transformation(extent = {{-86, -10}, {-66, 10}})));
        Interfaces.Flange_a flange_a1 annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
        Interfaces.FlangeChangeRight flangeChangeRight annotation(Placement(transformation(extent = {{66, -10}, {86, 10}})));
        Interfaces.Flange_b flange_b1 annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
        NominalLifeRating nominalLifeRating(p = p, C_100 = C_100) annotation(Placement(transformation(extent = {{20, -60}, {40, -40}})));
        Modelica.Blocks.Interfaces.RealOutput L_h
          "Connector of Real output signal"                                         annotation(Placement(transformation(extent = {{100, -60}, {120, -40}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {60, -110})));
        Modelica.Blocks.Interfaces.RealOutput P_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, -90}, {120, -70}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-60, -110})));
        HelpBlocks.AbsMax absMax
          annotation (Placement(transformation(extent={{-24,70},{-4,90}})));
        HelpBlocks.AbsMax absMax1
          annotation (Placement(transformation(extent={{-24,30},{-4,50}})));
        HelpBlocks.AbsMax absMax2
          annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
      equation
        connect(flangeChangeLeft.flange_a, flange_a1) annotation(Line(points = {{-81, 0}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(flangeChangeRight.flange_b, flange_b1) annotation(Line(points = {{81, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(flangeChangeLeft.flange_ax, flangeChangeRight.flange_ax) annotation(Line(points = {{-71, 10}, {71, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(accSensor.flange, speedSensor.flange) annotation(Line(points = {{-60, 80}, {-60, 40}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flangeChangeRight.flange_ax) annotation(Line(points = {{-60, 40}, {-60, 10}, {71, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(flangeChangeLeft.flange_N, forceSensor.flange_a) annotation(Line(points = {{-71, -10}, {-8, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.flange_b, flangeChangeRight.flange_N) annotation(Line(points = {{12, -10}, {71, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.f, nominalLifeRating.u1) annotation(Line(points = {{-6, -21}, {-6, -44}, {18, -44}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.v, nominalLifeRating.u2) annotation(Line(points = {{-39, 40}, {-39, -56}, {18, -56}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(nominalLifeRating.y, L_h) annotation(Line(points = {{41, -50}, {110, -50}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(accSensor.a, absMax.u) annotation (Line(
            points={{-39,80},{-26,80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax.y, a_max) annotation (Line(
            points={{-3,80},{110,80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(speedSensor.v, absMax1.u) annotation (Line(
            points={{-39,40},{-26,40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax1.y, v_max) annotation (Line(
            points={{-3,40},{110,40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax2.y, P_max) annotation (Line(
            points={{41,-80},{110,-80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMax2.u, nominalLifeRating.u1) annotation (Line(
            points={{18,-80},{2,-80},{2,-44},{18,-44}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-74, -60}, {66, 20}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Polygon(points=  {{-4, -40}, {-14, -16}, {6, -16}, {-4, -40}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-4, 0}, {-4, -16}}, color=  {0, 0, 0}), Line(points=  {{-74, 0}, {-4, 0}}, color=  {0, 0, 0}), Line(points=  {{-54, -40}, {-54, -60}}, color=  {0, 0, 0}), Line(points=  {{-34, -40}, {-34, -60}}, color=  {0, 0, 0}), Line(points=  {{-14, -40}, {-14, -60}}, color=  {0, 0, 0}), Line(points=  {{6, -40}, {6, -60}}, color=  {0, 0, 0}), Line(points=  {{26, -40}, {26, -60}}, color=  {0, 0, 0}), Line(points=  {{46, -40}, {46, -60}}, color=  {0, 0, 0}), Line(points=  {{-90, -80}, {-10, -80}}, color=  {0, 0, 0}), Polygon(points=  {{16, -82}, {-14, -72}, {-14, -92}, {16, -82}}, lineColor=  {128, 128, 128}, fillColor=  {128, 128, 128},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-74, 0}, {-94, 0}}, color=  {0, 0, 0}), Line(points=  {{66, 0}, {96, 0}}, color=  {0, 0, 127}), Text(extent=  {{-154, 100}, {146, 60}}, textString=  "%name", lineColor=  {0, 0, 255}), Line(points=  {{-74, 0}, {-94, 0}}, color=  {0, 0, 0})}));
      end LinearRollerGuideSensor;

      model NominalLifeRating
        "Calculation of nominal life rating for linear roller guides"
        extends Modelica.Blocks.Interfaces.SI2SO;
        parameter Real p "Life expectancy exponent";
        SI.Force P "Actual Equivalent load";
        SI.Force P_m "Equivalent load of cycle";
        parameter SI.Force C_100
          "Dynmamic load rating for a nominal lifetime 100 km";
        SI.AngularVelocity v "Velocity of slide";
        SI.AngularVelocity v_m "Average velocity";
        Real temp_P(start = 0);
        Real temp_v(start = 0);
        SI.Distance L_10 "Nominal life expectancy";
        SI.Conversions.NonSIunits.Time_hour L_h
          "Nominal life expectancy in hours";
      initial equation
        L_10 = 0;
        L_h = 0;
      equation
        P = abs(u1);
        v = abs(u2);
        der(temp_P) = P ^ p * v;
        der(temp_v) = v;
        when terminal() then
          P_m = (temp_P / temp_v) ^ (1 / p);
          v_m = temp_v / time;
          L_10 = (C_100 / P_m) ^ p * 10 ^ 5;
          L_h = L_10 / (60 * v_m);
        end when;
        y = L_h;
        annotation(Documentation(info = "
<HTML>
<p>
This blocks computes output <b>y</b> as <i>sum</i> of the
two input signals <b>u1</b> and <b>u2</b>:
</p>
<pre>
    <b>y</b> = k1*<b>u1</b> + k2*<b>u2</b>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   k1= +2, k2= -3

  results in the following equations:

     y = 2 * u1 - 3 * u2
</pre>

</HTML>
"),       Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Text(
                      extent={{-98,-52},{7,-92}},
                      lineColor={0,0,0},
                      textString="%k2"),Text(
                      extent={{-100,90},{5,50}},
                      lineColor={0,0,0},
                      textString="%k1"),Text(
                      extent={{-150,150},{150,110}},
                      textString="%name",
                      lineColor={0,0,255}),Line(points={{-100,60},{-40,60},{-30,
                40}}, color={0,0,255}),Ellipse(extent={{-50,50},{50,-50}},
                lineColor={0,0,255}),Line(points={{-100,-60},{-40,-60},{-30,-40}},
                color={0,0,255}),Line(points={{-15,-25.99},{15,25.99}}, color={
                0,0,0}),Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,127},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Line(points={{50,0},{100,0}},
                color={0,0,255}),Line(points={{-100,60},{-74,24},{-44,24}},
                color={0,0,127}),Line(points={{-100,-60},{-74,-28},{-42,-28}},
                color={0,0,127}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,127}),Line(points={{50,0},{100,0}}, color={0,0,127}),Text(
                      extent={{-38,34},{38,-34}},
                      lineColor={0,0,0},
                      textString="L_h"),Text(
                      extent={{-100,52},{5,92}},
                      lineColor={0,0,0},
                      textString="P"),Text(
                      extent={{-100,-52},{5,-92}},
                      lineColor={0,0,0},
                      textString="v")}),
          Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,255},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Text(
                      extent={{-98,-52},{7,-92}},
                      lineColor={0,0,0},
                      textString="%k2"),Text(
                      extent={{-100,90},{5,50}},
                      lineColor={0,0,0},
                      textString="%k1"),Line(points={{-100,60},{-40,60},{-30,40}},
                color={0,0,255}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,255}),Line(points={{-100,-60},{-40,-60},{-30,-40}}, color=
                {0,0,255}),Line(points={{-15,-25.99},{15,25.99}}, color={0,0,0}),
                Rectangle(
                      extent={{-100,-100},{100,100}},
                      lineColor={0,0,127},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Line(points={{50,0},{100,0}},
                color={0,0,255}),Line(points={{-100,60},{-74,24},{-44,24}},
                color={0,0,127}),Line(points={{-100,-60},{-74,-28},{-42,-28}},
                color={0,0,127}),Ellipse(extent={{-50,50},{50,-50}}, lineColor=
                {0,0,127}),Line(points={{50,0},{100,0}}, color={0,0,127}),Text(
                      extent={{-100,52},{5,92}},
                      lineColor={0,0,0},
                      textString="p"),Text(
                      extent={{-100,-52},{5,-92}},
                      lineColor={0,0,0},
                      textString="v"),Text(
                      extent={{-52,-20},{53,20}},
                      lineColor={0,0,0},
                      textString="L_h")}));
      end NominalLifeRating;

      type RailSystemType = enumeration(
          BallRailSystem "Life expectancy exponent p = 3",
          RollerRailSystem "Life expectancy exponent p = 10/3")                                                                                annotation(Evaluate = true);
    end Sensors;
  end TranslationalComponents;

  package RotationalComponents
    package BasicComponents


      model EfficiencyFactor
        "Calculates the friction torque by using the efficiency factor. Atan function is used to obtain a continous model. Static friction is zero"

        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
        parameter Real eta = 1 "efficiency factor";
        SI.AngularVelocity w "Angular velocity of flange";
        SI.Torque T_F(start = 0) "Friction torque";
        parameter SI.AngularVelocity w_ref = 10
          "Angular velocity after which the friction torque can is approx. equal to (1-eta)*flange_a.tau";
      equation
        flange_a.phi = flange_b.phi;
        w = der(flange_a.phi);
        T_F = -2 / Constants.pi * (1 - eta) * sign(flange_a.tau) * flange_a.tau * atan(w / w_ref);
        flange_a.tau + flange_b.tau + T_F = 0;
        annotation(Icon(graphics={  Rectangle(extent = {{-100, 10}, {54, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                      FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Rectangle(extent = {{50, 10}, {100, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Line(points = {{-80, -25}, {-60, -25}}, color = {0, 0, 0}), Line(points = {{60, -25}, {80, -25}}, color = {0, 0, 0}), Line(points = {{-70, -25}, {-70, -70}}, color = {0, 0, 0}), Line(points = {{70, -25}, {70, -70}}, color = {0, 0, 0}), Line(points = {{-80, 25}, {-60, 25}}, color = {0, 0, 0}), Line(points = {{60, 25}, {80, 25}}, color = {0, 0, 0}), Line(points = {{-70, 45}, {-70, 25}}, color = {0, 0, 0}), Line(points = {{70, 45}, {70, 25}}, color = {0, 0, 0}), Line(points = {{-70, -70}, {70, -70}}, color = {0, 0, 0}), Text(extent = {{-160, 140}, {140, 100}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-150, -80}, {150, -120}}, lineColor = {0, 0, 0}, textString = "eta=%eta"), Polygon(points = {{-99, 50}, {-70, 50}, {-70, 68}, {-80, 68}, {-60, 84}, {-40, 68}, {-50, 68}, {-50, 30}, {-99, 30}, {-99, 50}}, lineColor = {0, 0, 0}, fillColor = {255, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid)}));
      end EfficiencyFactor;

      model Bearing
        "Model of angular velocity dependent friction losses. Atan function is used to obtain a continous model. Static friction is zero "
        extends Modelica.Electrical.Machines.Interfaces.FlangeSupport;
        extends Modelica.Mechanics.Translational.Interfaces.PartialRigid(L = 0, s(start = 0));
        constant Real eps = 0.1;
        parameter SI.CoefficientOfFriction mu = 0.006 "Friction coefficient";
        parameter SI.Diameter d "Diameter of Bearing";
        SI.Force F_ax "Axial force on bearing";
        parameter SI.Force F_pr = 0 "Preload force";
        parameter SI.AngularVelocity w_ref = 10
          "Angular velocity after which the friction torque can is approx. equal to (1-eta)*flange_a.tau";
      equation
        tau = -2 / Constants.pi * abs(F_ax) * mu * (d / 2) * atan(w / w_ref);
        //w = der(flange.phi);
        F_ax = flange_a.f + F_pr;
        //forceFlange.f;
        0 = flange_b.f + flange_a.f;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-60, 60}, {60, -60}}, lineColor = {0, 0, 0}, fillColor = {175, 175, 175},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{-50, 50}, {50, -50}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{-12, 50}, {8, 30}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{-10, -30}, {10, -50}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{24, -10}, {44, -30}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{22, 34}, {42, 14}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{-44, 30}, {-24, 10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{-44, -12}, {-24, -32}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {135, 135, 135}), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {0, 0, 0}, fillColor = {175, 175, 175},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{-20, 20}, {20, -20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{-260, 101}, {40, 61}}, lineColor = {0, 0, 255}, textString = "spindle
bearing"), Rectangle(extent = {{-120, 20}, {-80, -20}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                  fillPattern =                                                                                                 FillPattern.Solid)}), Documentation(info = "<html>
<p>
The friction losses are considered by the equations
</p>
<pre>
   tau / tauRef = (+w / wRef) ^ power_w    for w &gt; +wLinear
 - tau / tauRef = (-w / wRef) ^ power_w    for w &lt; -wLinear
</pre>
<p>
with
</p>
<pre>
  tauRef * wRef = PRef
</pre>
<p>
being the friction torque at the referenc angular velocity
<code>wRef</code>. The exponent <code>power_w</code> is
approximately 1.5 for axial ventilation and approximately 2.0 for radial ventilation.
</p>
<p>
For stability reasons the friction torque <code>tau</code> is approximated by a linear curve
</p>
<pre>
  tau / tauLinear = w / wLinear
</pre>
<p>
with
</p>
<pre>
  tauLinear = tauRef*(wLinear/wRef) ^ power_w
</pre>
<p>
in the range <code> -wLinear &le; w &le; wLinear</code> with <code>wLinear = 0.001 * wRef</code>. The relationship of torque
and angular velocity is depicted in Fig. 1
</p>
<table border=0 cellspacing=0 cellpadding=1>
  <tr><td> <img src=\"modelica://Modelica/Resources/Images/Electrical/Machines/frictiontorque.png\"> </td>
  </tr>
  <tr><td> <b> Fig. 1: </b>Friction loss torque versus angular velocity for <code>power_w = 2</code></td>
  </tr>
</table>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Electrical.Machines.Losses.FrictionParameters\">FrictionParameters</a>
</p>
<p>
If it is desired to neglect friction losses, set <code>frictionParameters.PRef = 0</code> (this is the default).
</p>
</html>"));
      end Bearing;

      model Torque_PT1
        parameter SI.Time T = 0.01 "Time for PT1";
        Modelica.Mechanics.Rotational.Sources.Torque torque annotation(Placement(transformation(extent = {{20, -10}, {40, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange1
          "Flange of shaft"                                                         annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Continuous.FirstOrder firstOrder(k = 1, T = T) annotation(Placement(transformation(extent = {{-40, -10}, {-20, 10}})));
        Modelica.Blocks.Interfaces.RealInput tau
          "Accelerating torque acting at flange (= -flange.tau)"                                        annotation(Placement(transformation(extent = {{-140, -20}, {-100, 20}}, rotation = 0), iconTransformation(extent = {{-130, -20}, {-90, 20}})));
      equation
        connect(firstOrder.y, torque.tau) annotation(Line(points = {{-19, 0}, {18, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(firstOrder.u, tau) annotation(Line(points = {{-42, 0}, {-120, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(flange1, torque.flange) annotation(Line(points = {{100, 0}, {40, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(visible = not useSupport, points = {{-40, -120}, {-20, -100}}, color = {0, 0, 0}), Line(visible = not useSupport, points = {{-20, -120}, {0, -100}}, color = {0, 0, 0}), Line(visible = not useSupport, points = {{0, -120}, {20, -100}}, color = {0, 0, 0}), Line(visible = not useSupport, points = {{20, -120}, {40, -100}}, color = {0, 0, 0}), Line(visible = not useSupport, points = {{-20, -100}, {40, -100}}, color = {0, 0, 0}), Text(extent = {{-140, 110}, {160, 70}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-52, -29}, {-131, -70}}, lineColor = {0, 0, 0}, textString = "tau"), Line(points = {{-78, 0}, {-54, 30}, {-26, 52}, {8, 62}, {38, 56}, {58, 44}, {74, 28}, {86, 14}, {96, 0}}, color = {0, 0, 0}, thickness = 0.5), Polygon(points = {{96, 0}, {76, 58}, {47, 27}, {96, 0}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-20, -30}, {40, -30}}, color = {0, 0, 0}), Line(points = {{10, -30}, {10, -101}}, color = {0, 0, 0}), Line(points = {{-20, -50}, {0, -30}}, color = {0, 0, 0}), Line(points = {{0, -50}, {20, -30}}, color = {0, 0, 0}), Line(points = {{20, -50}, {40, -30}}, color = {0, 0, 0}), Line(points = {{-44, -42}, {-28, -28}, {-6, -16}, {14, -14}, {32, -18}, {46, -26}, {58, -36}, {66, -46}, {74, -58}}, color = {0, 0, 0}), Polygon(points = {{-51, -66}, {-34, -42}, {-48, -36}, {-51, -66}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid)}));
      end Torque_PT1;


      model NoLoad_Friction
        "No load friction is independent of the load and a standard parameter for ball screw drives"
        SI.AngularVelocity w "Velocity of guide";
        parameter SI.Torque M_Fr "Friction force";
        SI.Torque M_Fr_dir "Directed friction force";
        parameter SI.AngularVelocity w_ref = 300
          "Velocity after which the friction force is approx. equal to mu*F_N";
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
      equation
        flange_a.phi = flange_b.phi;
        w = der(flange_a.phi);
        M_Fr_dir = -2 / Constants.pi * M_Fr * atan(w / w_ref);
        flange_a.tau + flange_b.tau + M_Fr = 0;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 10}, {54, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Rectangle(extent = {{50, 10}, {100, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Line(points = {{-80, -25}, {-60, -25}}, color = {0, 0, 0}), Line(points = {{60, -25}, {80, -25}}, color = {0, 0, 0}), Line(points = {{-70, -25}, {-70, -70}}, color = {0, 0, 0}), Line(points = {{70, -25}, {70, -70}}, color = {0, 0, 0}), Line(points = {{-80, 25}, {-60, 25}}, color = {0, 0, 0}), Line(points = {{60, 25}, {80, 25}}, color = {0, 0, 0}), Line(points = {{-70, 45}, {-70, 25}}, color = {0, 0, 0}), Line(points = {{70, 45}, {70, 25}}, color = {0, 0, 0}), Line(points = {{-70, -70}, {70, -70}}, color = {0, 0, 0}), Text(extent = {{-150, -82}, {150, -122}}, lineColor = {0, 0, 0}, textString = "M_Fr=%M_Fr"), Polygon(points = {{-99, 50}, {-70, 50}, {-70, 68}, {-80, 68}, {-60, 84}, {-40, 68}, {-50, 68}, {-50, 30}, {-99, 30}, {-99, 50}}, lineColor = {0, 0, 0}, fillColor = {255, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{-220, 140}, {220, 100}}, lineColor = {0, 0, 255}, textString = "%name")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
      end NoLoad_Friction;

    end BasicComponents;

    package TimingBelt
      model BeltPulley "Model for a toothed belt drive pulley"
        parameter SI.Diameter d_eff(displayUnit = "mm") = 0.05602
          "Effective diameter"                                                         annotation(Dialog(group = "Physical constants"));
        parameter Integer z = 22 "Number of teeth" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inertia J = (d_eff / 0.08149) ^ 4 * 0.000655
          "Moment of inertia of pulley"                                                         annotation(Dialog(group = "Physical constants"));
        extends
          Modelica.Mechanics.Rotational.Interfaces.PartialElementaryRotationalToTranslational;
      protected
        parameter StateSelect stateSelect = StateSelect.default
          "Priority to use phi and w as states"                                                       annotation(HideResult = true, Dialog(tab = "Advanced"));
      public
        SI.Angle phi(stateSelect = stateSelect)
          "Absolute rotation angle of component"                                       annotation(Dialog(group = "Initialization", showStartAttribute = true));
        SI.AngularVelocity w(stateSelect = stateSelect)
          "Absolute angular velocity of component (= der(phi))"                                               annotation(Dialog(group = "Initialization", showStartAttribute = true));
        SI.AngularAcceleration a
          "Absolute angular acceleration of component (= der(w))";
      equation
        (flangeR.phi - internalSupportR.phi) * (d_eff / 2) = flangeT.s - internalSupportT.s;
        phi = flangeR.phi;
        w = der(phi);
        a = der(w);
        J * a = d_eff / 2 * flangeT.f + flangeR.tau;
        annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(visible = not useSupportT, points = {{85, -110}, {95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{95, -110}, {105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{105, -110}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{85, -100}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -110}, {-105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-105, -110}, {-95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-95, -110}, {-85, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -100}, {-85, -100}}, color = {0, 0, 0}), Rectangle(extent = {{-98, 10}, {-44, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Ellipse(extent = {{-48, 80}, {12, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {160, 160, 164}), Rectangle(extent = {{-20, 80}, {12, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Ellipse(extent = {{-16, 80}, {44, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {160, 160, 164}), Ellipse(extent = {{-2, 52}, {34, -52}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{12, 10}, {20, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {192, 192, 192}), Rectangle(extent = {{16, 10}, {60, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Ellipse(extent = {{56, 10}, {64, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {192, 192, 192}), Text(extent = {{-150, 140}, {150, 100}}, textString = "%name", lineColor = {0, 0, 255}), Polygon(points = {{80, 10}, {80, 26}, {60, 26}, {60, 20}, {70, 20}, {70, -20}, {60, -20}, {60, -26}, {80, -26}, {80, -10}, {90, -10}, {90, 10}, {80, 10}}, lineColor = {0, 127, 0}, fillColor = {0, 127, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-100, -20}, {-60, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-100, -20}, {-100, -100}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-100, 20}, {-60, 20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{100, -90}, {-40, -90}}, color = {0, 127, 0}, smooth = Smooth.None), Line(points = {{70, -26}, {70, -50}, {100, -50}, {100, -100}}, color = {0, 127, 0}, smooth = Smooth.None)}));
      end BeltPulley;

      model PulleyConstraints
        "Model for a toothed belt drive pully with constraints in force and power"
        parameter SI.Diameter d_eff(displayUnit = "mm") = 0.05602
          "Effective diameter"                                                           annotation(Dialog(group = "Physical constants"));
        parameter Integer z = 22 "Number of teeth" annotation(Dialog(group = "Physical constants"));
        parameter SI.Inertia J = 1.44e-4 "Moment of inertia of pulley" annotation(Dialog(group = "Physical constants"));
        parameter Real eta = 0.98 "Efficiency factor" annotation(Dialog(group = "Physical constants"));
        parameter Real exponent = 0.80919
          "Exponent for permissable power table (P_perm[W] = factor*n[U/min]^exponent)"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter Real factor = 13.17
          "Factor for perissable power table (P_perm[W] = factor*n[U/min]^exponent)"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_perm = 2012 "Permissable force" annotation(Dialog(group = "Limiting data"));
        parameter Real S_Operation = 1.3 "See belt catalogue";
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flangeR
          "Flange of rotational shaft"                                                         annotation(Placement(transformation(extent={{-110,
                  -10},{-90,10}},                                                                                                    rotation = 0), iconTransformation(extent={{-110,
                  -10},{-90,10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flangeT
          "Flange of translational rod"                                                            annotation(Placement(transformation(extent={{90,10},
                  {110,-10}},                                                                                                    rotation = 0), iconTransformation(extent={{90,10},
                  {110,-10}})));
        BeltPulley beltPulley(d_eff = d_eff, z = z, J = J) annotation(Placement(transformation(extent = {{-10, 18}, {10, 38}})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor annotation(Placement(transformation(extent = {{-80, 18}, {-60, 38}})));
        Modelica.Blocks.Math.Abs abs1 annotation(Placement(transformation(extent = {{-50, -50}, {-30, -30}})));
        Modelica.Blocks.Math.Gain gain(k = S_Operation) annotation(Placement(transformation(extent = {{-6, -50}, {14, -30}})));
        Modelica.Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{42, -50}, {62, -30}})));
        HelpBlocks.PolynomialFunction perm_power(exponent = exponent, factor = factor) annotation(Placement(transformation(extent = {{0, -90}, {20, -70}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-58, -90}, {-38, -70}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm speed_to_rpm annotation(Placement(transformation(extent = {{-30, -90}, {-10, -70}})));
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{28, 18}, {48, 38}})));
        Real con_F "Constraint for force";
        Real con_P "Contraitnt for power";
        Real util_F "Utilization for force";
        Real util_P "Utilization for power";
        BasicComponents.EfficiencyFactor efficiencyFactor(eta = eta, w_ref = 10) annotation(Placement(transformation(extent = {{-46, 18}, {-26, 38}})));
        HelpBlocks.Max max_con_P
          annotation (Placement(transformation(extent={{78,-50},{98,-30}})));
        HelpBlocks.Max max_perm
          annotation (Placement(transformation(extent={{78,-90},{98,-70}})));
        HelpBlocks.AbsMax absMax_Force
          annotation (Placement(transformation(extent={{72,70},{92,90}})));
      initial equation
        util_F = 0;
        util_P = 0;
      algorithm
        when terminal() then
          util_F := absMax_Force.y / F_perm;
          util_P := 1 + max_con_P.y / max_perm.y;
        end when;
      equation
        con_F = absMax_Force.y - F_perm;
        con_P = max_con_P.y;
        connect(powerSensor.power, abs1.u) annotation(Line(points = {{-78, 17}, {-78, -40}, {-52, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(abs1.y, gain.u) annotation(Line(points = {{-29, -40}, {-8, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(gain.y, feedback.u1) annotation(Line(points = {{15, -40}, {44, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(flangeR, powerSensor.flange_a) annotation(Line(points={{-100,0},
                {-90,0},{-90,28},{-80,28}},                                                      color = {0, 0, 0}, smooth = Smooth.None));
        connect(flangeR, speedSensor.flange) annotation(Line(points={{-100,0},{
                -80,0},{-80,-80},{-58,-80}},                                                                           color = {0, 0, 0}, smooth = Smooth.None));
        connect(speedSensor.w, speed_to_rpm.u) annotation(Line(points = {{-37, -80}, {-32, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speed_to_rpm.y, perm_power.u) annotation(Line(points = {{-9, -80}, {-2, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(beltPulley.flangeT, forceSensor.flange_a) annotation(Line(points = {{10, 28}, {28, 28}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.flange_b, flangeT) annotation(Line(points={{48,28},
                {74,28},{74,0},{100,0}},                                                       color = {0, 127, 0}, smooth = Smooth.None));
        connect(powerSensor.flange_b, efficiencyFactor.flange_a) annotation(Line(points = {{-60, 28}, {-46, 28}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(efficiencyFactor.flange_b, beltPulley.flangeR) annotation(Line(points = {{-26, 28}, {-10, 28}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(feedback.y, max_con_P.u) annotation (Line(
            points={{61,-40},{76,-40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(feedback.u2, perm_power.y) annotation (Line(
            points={{52,-48},{52,-80},{21,-80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(max_perm.u, perm_power.y) annotation (Line(
            points={{76,-80},{21,-80}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(forceSensor.f, absMax_Force.u) annotation (Line(
            points={{30,17},{50,17},{50,80},{70,80}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Icon(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                        graphics={  Line(visible = not useSupportT, points = {{85, -110}, {95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{95, -110}, {105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{105, -110}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{85, -100}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -110}, {-105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-105, -110}, {-95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-95, -110}, {-85, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -100}, {-85, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{85, -110}, {95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{95, -110}, {105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{105, -110}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportT, points = {{85, -100}, {115, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -110}, {-105, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-105, -110}, {-95, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-95, -110}, {-85, -100}}, color = {0, 0, 0}), Line(visible = not useSupportR, points = {{-115, -100}, {-85, -100}}, color = {0, 0, 0}), Rectangle(extent = {{-98, 10}, {-44, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Ellipse(extent = {{-48, 80}, {12, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {160, 160, 164}), Rectangle(extent = {{-20, 80}, {12, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Ellipse(extent = {{-16, 80}, {44, -80}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {160, 160, 164}), Ellipse(extent = {{-2, 52}, {34, -52}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                  fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{12, 10}, {20, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {192, 192, 192}), Rectangle(extent = {{16, 10}, {60, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Ellipse(extent = {{56, 10}, {64, -10}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                   FillPattern.Sphere, fillColor = {192, 192, 192}), Text(extent = {{-150, 140}, {150, 100}}, textString = "%name", lineColor = {0, 0, 255}), Polygon(points = {{80, 10}, {80, 26}, {60, 26}, {60, 20}, {70, 20}, {70, -20}, {60, -20}, {60, -26}, {80, -26}, {80, -10}, {90, -10}, {90, 10}, {80, 10}}, lineColor = {0, 127, 0}, fillColor = {0, 127, 0},
                  fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-100, -20}, {-60, -20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-100, -20}, {-100, -100}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{-100, 20}, {-60, 20}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{70, -26}, {70, -50}, {100, -50}, {100, -100}}, color = {0, 127, 0}, smooth = Smooth.None)}), Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                                                    graphics));
      end PulleyConstraints;

      model TimingBeltDrive
        "Model of a toothed belt drive. Can be used with data from companies such as Walther Flender"
        parameter Boolean elastic = false "True, if timing belt is elastic";
        parameter SI.Diameter d_eff_s(displayUnit = "mm") "Effective diameter" annotation (Dialog(group = "Small pulley"));
        parameter Integer z_s "Number of teeth" annotation(Dialog(group = "Small pulley"));
        parameter SI.Inertia J_s "Moment of inertia of pulley" annotation(Dialog(group = "Small pulley"));
        parameter SI.Diameter d_eff_l(displayUnit = "mm") "Effective diameter" annotation (Dialog(group = "Large pulley"));
        parameter Integer z_l "Number of teeth" annotation(Dialog(group = "Large pulley"));
        parameter SI.Inertia J_l "Moment of inertia of pulley" annotation(Dialog(group = "Large pulley"));
        parameter SI.Length width(displayUnit = "mm") "Width of belt" annotation(Dialog(group = "Belt"));
        parameter SI.Length standard_width(displayUnit = "mm")
          "Width at which permissable power is given"                                  annotation(Dialog(group = "Belt"));
        parameter Real exponent
          "Exponent for permissable power table (P_perm[W] = factor*n[U/min]^exponent)"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter Real factor
          "Factor for perissable power table (P_perm[W] = factor*n[U/min]^exponent)"
                                                                                                         annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_perm "Permissable force" annotation(Dialog(group = "Limiting data"));
        parameter SI.TranslationalSpringConstant c(final min = 0) = 1e6
          "Spring constant ";
        parameter Real S_Operation = 1.3 "See belt catalogue";
        parameter Real eta = 0.95 "Efficiency factor";
        BeltPulley largePulley(d_eff = d_eff_l, z = z_l, J = J_l) annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {70, 0})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flangeR1
          "Flange of rotational shaft"                                                          annotation(Placement(transformation(extent = {{-130, -10}, {-110, 10}}), iconTransformation(extent = {{-130, -10}, {-110, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flangeR2
          "Flange of rotational shaft"                                                          annotation(Placement(transformation(extent = {{110, -10}, {130, 10}}), iconTransformation(extent = {{110, -10}, {130, 10}})));
        BeltPulley smallPulley(d_eff = d_eff_s, z = z_s, J = J_s) annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor annotation(Placement(transformation(extent = {{-100, -10}, {-80, 10}})));
        Modelica.Blocks.Math.Abs abs1 annotation(Placement(transformation(extent = {{-80, -50}, {-60, -30}})));
        Modelica.Blocks.Math.Gain g_Operation(k = S_Operation) annotation(Placement(transformation(extent = {{-26, -50}, {-6, -30}})));
        Modelica.Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{22, -70}, {42, -50}})));
        Modelica.Blocks.Math.Gain g_width(k = width / standard_width) annotation(Placement(transformation(extent = {{10, -90}, {30, -70}})));
        HelpBlocks.PolynomialFunction perm_power(exponent = exponent, factor = factor) annotation(Placement(transformation(extent = {{-20, -90}, {0, -70}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-80, -90}, {-60, -70}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm speed_to_rpm annotation(Placement(transformation(extent = {{-50, -90}, {-30, -70}})));
        BasicComponents.EfficiencyFactor efficiencyFactor(eta = eta, w_ref = 10) annotation(Placement(transformation(extent = {{-66, -10}, {-46, 10}})));
        TimingBelt_elastic timingBelt_elastic(c = c) annotation(Placement(transformation(extent = {{20, -10}, {40, 10}})));
        HelpBlocks.Max max_con_P
          annotation (Placement(transformation(extent={{60,-70},{80,-50}})));
        HelpBlocks.AbsMax absMax_Force
          annotation (Placement(transformation(extent={{60,-42},{80,-22}})));
      equation
        connect(powerSensor.power, abs1.u) annotation(Line(points = {{-98, -11}, {-98, -40}, {-82, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(abs1.y, g_Operation.u) annotation(Line(points = {{-59, -40}, {-28, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(g_Operation.y, feedback.u1) annotation(Line(points = {{-5, -40}, {10, -40}, {10, -60}, {24, -60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(g_width.y, feedback.u2) annotation(Line(points = {{31, -80}, {32, -80}, {32, -68}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, speed_to_rpm.u) annotation(Line(points = {{-59, -80}, {-52, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speed_to_rpm.y, perm_power.u) annotation(Line(points = {{-29, -80}, {-22, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(perm_power.y, g_width.u) annotation(Line(points = {{1, -80}, {8, -80}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(powerSensor.flange_b, efficiencyFactor.flange_a) annotation(Line(points = {{-80, 0}, {-66, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(efficiencyFactor.flange_b, smallPulley.flangeR) annotation(Line(points = {{-46, 0}, {-20, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(flangeR1, flangeR1) annotation(Line(points = {{-120, 0}, {-120, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(largePulley.flangeR, flangeR2) annotation(Line(points = {{80, 0}, {120, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(powerSensor.flange_a, flangeR1) annotation(Line(points = {{-100, 0}, {-112, 0}, {-112, 0}, {-120, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flangeR1) annotation(Line(points = {{-80, -80}, {-108, -80}, {-108, 0}, {-120, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(smallPulley.flangeT, timingBelt_elastic.flange_a1) annotation(Line(points = {{0, 0}, {20, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(largePulley.flangeT, timingBelt_elastic.flange_b1) annotation(Line(points = {{60, 1.22125e-015}, {50, 1.22125e-015}, {50, 0}, {40, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(feedback.y, max_con_P.u) annotation (Line(
            points={{41,-60},{58,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(timingBelt_elastic.f1, absMax_Force.u) annotation (Line(
            points={{31.2,-4.6},{31.2,-31.3},{58,-31.3},{58,-32}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-150,
                  -100},{150,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-150, -100}, {150, 100}}), graphics={  Ellipse(extent=  {{-80, 20}, {-40, -20}}, lineColor=  {0, 0, 0}), Ellipse(extent=  {{0, 40}, {80, -40}}, lineColor=  {0, 0, 0}), Ellipse(extent=  {{-62, 2}, {-58, -2}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{38, 2}, {42, -2}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-60, 20}, {34, 40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-60, -20}, {34, -40}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{-80, 0}, {-118, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{110, 0}, {80, 0}}, color=  {0, 0, 0}, smooth=  Smooth.None), Text(extent=  {{-150, 100}, {150, 60}}, lineColor=  {0, 0, 255}, textString=  "%name"), Text(extent=  {{-150, -40}, {150, -80}}, lineColor=  {0, 0, 255}, textString=  "z1 = %z_s"), Text(extent=  {{-150, -74}, {150, -114}}, lineColor=  {0, 0, 255}, textString=  "z2 = %z_l")}));
      end TimingBeltDrive;

      model TimingBelt_elastic
        parameter SI.TranslationalSpringConstant c(final min = 0) = 1e6
          "Spring constant ";
        parameter Boolean elastic = false "True, if timing belt is elastic";
        Modelica.Mechanics.Translational.Components.Spring spring(c = c) if elastic annotation(Placement(transformation(extent = {{-30, -10}, {-10, 10}})));
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor annotation(Placement(transformation(extent = {{10, -10}, {30, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a1
          "Left flange of compliant 1-dim. translational component"                                                              annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b1
          "(right) driven flange (flange axis directed out of cut plane)"                                                              annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Interfaces.RealOutput f1
          "Force in flange_a and flange_b (f = flange_a.f = -flange_b.f)"                                        annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {12, -46})));
      equation
        if elastic then
          connect(spring.flange_a, flange_a1) annotation(Line(points = {{-30, 0}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
          connect(spring.flange_b, forceSensor.flange_a) annotation(Line(points = {{-10, 0}, {10, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        else
          connect(forceSensor.flange_a, flange_a1) annotation(Line(points = {{10, 0}, {10, 16}, {-100, 16}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        end if;
        connect(forceSensor.flange_b, flange_b1) annotation(Line(points = {{30, 0}, {100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(forceSensor.f, f1) annotation(Line(points = {{12, -11}, {12, -46}}, color = {0, 0, 127}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-70, 10}, {70, -10}}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{-70, -22}, {-50, -2}}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{-50, -20}, {-30, 0}}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{-30, -20}, {-10, 0}}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{-10, -20}, {10, 0}}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{10, -20}, {30, 0}}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{30, -20}, {50, 0}}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Ellipse(extent=  {{50, -20}, {70, 0}}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None), Line(points=  {{-96, 0}, {-70, 0}, {92, 0}}, pattern=  LinePattern.None, smooth=  Smooth.None), Text(extent=  {{-150, 100}, {150, 60}}, textString=  "%name", lineColor=  {0, 0, 255})}));
      end TimingBelt_elastic;
    end TimingBelt;

    package Sensors
      model SpeedSensor_rpm
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-58, -10}, {-38, 10}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm absMaxSpeed_to_rpm annotation(Placement(transformation(extent = {{20, -40}, {40, -20}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm meanSpeed_to_rpm annotation(Placement(transformation(extent = {{20, -70}, {40, -50}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm speed annotation(Placement(transformation(extent = {{-16, -10}, {4, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange
          "Flange of shaft from which sensor information shall be measured"                                                        annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Blocks.Interfaces.RealOutput actual
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{70, 16}, {90, 36}}), iconTransformation(extent = {{70, 16}, {90, 36}})));
        Modelica.Blocks.Interfaces.RealOutput max
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{70, -10}, {90, 10}}), iconTransformation(extent = {{70, -10}, {90, 10}})));
        Modelica.Blocks.Interfaces.RealOutput mean
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{70, -34}, {90, -14}}), iconTransformation(extent = {{70, -34}, {90, -14}})));
        HelpBlocks.AbsMax absMaxSpeed
          annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
        HelpBlocks.MeanValue meanSpeed
          annotation (Placement(transformation(extent={{-16,-70},{4,-50}})));
      equation
        connect(speedSensor.flange, flange) annotation(Line(points = {{-58, 0}, {-82, 0}, {-82, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(speed.y, actual) annotation(Line(points = {{5, 0}, {42, 0}, {42, 26}, {80, 26}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxSpeed_to_rpm.y, max) annotation(Line(points = {{41, -30}, {60, -30}, {60, 0}, {80, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(meanSpeed_to_rpm.y, mean) annotation(Line(points = {{41, -60}, {60, -60}, {60, -24}, {80, -24}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxSpeed.y, absMaxSpeed_to_rpm.u) annotation (Line(
            points={{5,-30},{18,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed_to_rpm.u, meanSpeed.y) annotation (Line(
            points={{18,-60},{5,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(speedSensor.w, speed.u) annotation (Line(
            points={{-37,0},{-18,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.u, speed.u) annotation (Line(
            points={{-18,-30},{-24,-30},{-24,0},{-18,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.u, speed.u) annotation (Line(
            points={{-18,-60},{-22,-60},{-22,-30},{-24,-30},{-24,0},{-18,0}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points=  {{-74, 0}, {-94, 0}}, color=  {0, 0, 0}), Text(extent=  {{146, 80}, {-154, 120}}, textString=  "%name", lineColor=  {0, 0, 255}), Ellipse(extent=  {{-74, 70}, {66, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-4, 70}, {-4, 40}}, color=  {0, 0, 0}), Line(points=  {{18.9, 32.8}, {36.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-26.9, 32.8}, {-44.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{33.6, 13.7}, {61.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-41.6, 13.7}, {-69.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-4, 0}, {5.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-4.48, 31.6}, {14, 26}, {14, 57.2}, {-4.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-9, 5}, {1, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{66, -30}, {116, -70}}, lineColor=  {0, 0, 0}, textString=  "rpm")}));
      end SpeedSensor_rpm;

      model TorqueSensor
        Modelica.Mechanics.Rotational.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Interfaces.RealOutput max
          "Connector of Real output signal"                                         annotation(Placement(transformation(extent = {{100, -40}, {120, -20}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-92, -110})));
        Modelica.Blocks.Interfaces.RealOutput rms
          "Connector of Real output signal"                                         annotation(Placement(transformation(extent = {{100, -80}, {120, -60}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-72, -110})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
        HelpBlocks.RootMeanSquareValue rootMeanSquareValue annotation(Placement(transformation(extent = {{34, -80}, {54, -60}})));
        HelpBlocks.AbsMax absMaxTorque
          annotation (Placement(transformation(extent={{34,-40},{54,-20}})));
      equation
        connect(torqueSensor.flange_a, flange_a) annotation(Line(points = {{-10, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.flange_b, flange_b) annotation(Line(points = {{10, 0}, {100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(rootMeanSquareValue.y, rms) annotation(Line(points = {{55, -70}, {110, -70}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxTorque.y, max) annotation (Line(
            points={{55,-30},{110,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.u, torqueSensor.tau) annotation (Line(
            points={{32,-30},{-8,-30},{-8,-11}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(rootMeanSquareValue.u, torqueSensor.tau) annotation (Line(
            points={{32,-70},{-8,-70},{-8,-11}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Text(extent=  {{-150, 73}, {150, 113}}, textString=  "%name", lineColor=  {0, 0, 255}), Ellipse(extent=  {{-70, 70}, {70, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 70}, {0, 40}}, color=  {0, 0, 0}), Line(points=  {{22.9, 32.8}, {40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-22.9, 32.8}, {-40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{37.6, 13.7}, {65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-37.6, 13.7}, {-65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{0, 0}, {9.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-0.48, 31.6}, {18, 26}, {18, 57.2}, {-0.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-5, 5}, {5, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-50, -80}, {50, -120}}, lineColor=  {0, 0, 0}, textString=  "tau"), Line(points=  {{-80, -100}, {-80, 0}}, color=  {0, 0, 127})}));
      end TorqueSensor;
    end Sensors;

    package Servogear
      model Servogear_SEW
        "Servogear according to SEW Eurodrive. PSF 221 EPH 02 used as reference"
        parameter Real ratio = 3
          "Transmission ratio (flange_a.phi/flange_b.phi)"                        annotation(Dialog(group = "Physical data"));
        parameter Real eta = 0.99 "Efficiency factor" annotation(Dialog(group = "Physical data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_e = 1500
          "Rated speed"                                                                  annotation(Dialog(group = "Engineering data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_e_pk = 7500
          "Maximum permitted input speed for short-time duty"                                                                     annotation(Dialog(group = "Limiting data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_a_k = 1890
          "Breakpoint speed (germ. Knickdrehzahl)"                                                                    annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_a_max = 40
          "Maximum permitted output torque for continous duty"                                annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_a_pk = 53
          "Maximum permitted output torque for short-time duty"                               annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_a_NOTAUS = 75
          "Maximum permitted output emergency stop torque, max 1000 stops"                                   annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_R_a_max = 2870
          "Maximum permitted overhung load at the output shaft for continous duty"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_R_a_pk = 2720
          "Maximum permitted overhung load at the output shaft for short-time duty"
                                                                                                              annotation(Dialog(group = "Limiting data"));
        parameter SI.Inertia J_GA = 0.73e-4
          "Moment of inertia of the gear unit including adapter with reference to the input shaft"
                                                                                                              annotation(Dialog(group = "Physical data"));
        parameter Real a_0 = 99.00
          "Gear unit constant 1 as regards to the rise in temperature in the gear unit"
                                                                                                              annotation(Dialog(group = "Engineering data"));
        parameter Real a_1 = -0.071
          "Gear unit constant 2 as regards to the rise in temperature in the gear unit"
                                                                                                              annotation(Dialog(group = "Engineering data"));
        parameter Real a_2 = 0
          "Gear unit constant 3 as regards to the rise in temperature in the gear unit"
                                                                                                             annotation(Dialog(group = "Engineering data"));
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = ratio) annotation(Placement(transformation(extent = {{-80, -10}, {-60, 10}})));
        BasicComponents.EfficiencyFactor efficiencyFactor(eta = eta, w_ref = 0.01 * n_e * 2 * Constants.pi / 60) annotation(Placement(transformation(extent = {{32, -10}, {52, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Flange of right shaft"                                                          annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-120, -10}, {-100, 10}}), iconTransformation(extent = {{-120, -10}, {-100, 10}})));
        Req_SEW s annotation(Placement(transformation(extent = {{-10, -10}, {10, 10}})));
        SI.Torque con_M_a_max
          "Maximum permitted output torque for continous duty";
        SI.Conversions.NonSIunits.AngularVelocity_rpm con_n_e_pk
          "Maximum permitted input speed for short-time duty";
        SI.Conversions.NonSIunits.AngularVelocity_rpm n_e_max
          "maximum speed of input shaft";
        Boolean dec_n_a_m
          "true, if mean speed is smaller than breakpoint speed (s.n_a_mean <= n_a_k)";
        SI.Torque temp_con_M_a_eff
          "Maximum permitted output torque for continous duty below breakpoint speed";
        SI.Torque temp_con_M_a_kub
          "Maximum permitted output torque for continous duty above breakpoint speed";
        Real f_k "Speed ratio";
        SI.Torque con_M_a "Maximum permitted output torque for continous duty";
        SI.Torque M_therm "Thermically admissible torque";
        SI.Torque con_M_therm "Thermically admissible torque";
      equation
        n_e_max = s.n_a_max * ratio;
        con_M_a_max = s.M_a_max - M_a_pk;
        con_n_e_pk = n_e_max - n_e_pk;
        dec_n_a_m = if s.n_a_mean <= n_a_k then true else false;
        temp_con_M_a_eff = s.M_a_8 - M_a_max;
        f_k = (s.n_a_mean / n_a_k) ^ 0.3;
        temp_con_M_a_kub = if dec_n_a_m then temp_con_M_a_eff else temp_con_M_a_kub;
        con_M_a = if dec_n_a_m then temp_con_M_a_eff else temp_con_M_a_kub;
        M_therm = if time <= 0 then a_0 else a_0 + a_1 * s.n_a_mean + a_2 / s.n_a_mean ^ 1.2;
        con_M_therm = s.M_a_1_2 - M_therm;
        connect(flange_b, efficiencyFactor.flange_b) annotation(Line(points = {{110, 0}, {52, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(flange_a, idealGear.flange_a) annotation(Line(points = {{-110, 0}, {-80, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(s.flange_b, efficiencyFactor.flange_a) annotation(Line(points = {{10, 0}, {32, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGear.flange_b, s.flange_a) annotation(Line(points = {{-60, 0}, {-10, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-100, 100}, {100, -100}}, fillColor=  {95, 95, 95},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None, lineColor=  {0, 0, 0}), Ellipse(extent=  {{-80, 80}, {80, -80}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-40, 40}, {40, -40}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-80, 20}, {-40, -20}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{40, 20}, {80, -20}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-53, -30}, {-13, -70}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{13, -30}, {53, -70}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-53, 70}, {-13, 30}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{13, 70}, {53, 30}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-84, 24}, {-36, -24}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-56, 74}, {-8, 26}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{9, 73}, {57, 25}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{36, 24}, {84, -24}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{9, -26}, {57, -74}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-57, -26}, {-9, -74}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-36, 35}, {36, -36}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-84, 83}, {84, -84}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash)}));
      end Servogear_SEW;

      model Servogear_Wittenstein "Servogear according to Wittenstein SP+"
        parameter Real ratio = 3
          "Transmission ratio (flange_a.phi/flange_b.phi)"                        annotation(Dialog(group = "Physical data"));
        parameter SI.Torque M_2_lim = 310
          "Max acceleration torque at output (max 1000 cycles per hour)"                                 annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_2_N_lim = 130
          "Maximum permitted output torque for short-time duty"                                   annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_a_emergency = 1000
          "Maximum permitted output emergency stop torque, max 1000 stops"                                        annotation(Dialog(group = "Limiting data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_1_N_lim = 3000
          "Nominal input speed"                                                                        annotation(Dialog(group = "Limiting data"));
        parameter SI.Conversions.NonSIunits.AngularVelocity_rpm n_1_lim = 6000
          "Maximum permitted input speed for short-time duty"                                                                      annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_noLoad_12 = 5.1 "mean no load running torque" annotation(Dialog(group = "Physical data"));
        parameter SI.RotationalSpringConstant c_t12(final min = 0, start = 1.0e5) = 53 * 3437.7
          "Torsional rigidity"                                                                                       annotation(Dialog(group = "Physical data"));
        parameter SI.Force F_2_a_max = 9870
          "Maximum axial force at output shaft"                                   annotation(Dialog(group = "Limiting data"));
        parameter SI.Force F_2_r_max = 9900
          "Maximum radial force at output shaft"                                   annotation(Dialog(group = "Limiting data"));
        parameter Real eta = 0.99 "Efficiency at full load" annotation(Dialog(group = "Physical data"));
        parameter SI.Inertia J_GA = 14.9e-4
          "Moment of inertia of the gear unit including adapter with reference to the input shaft"
                                                                                                              annotation(Dialog(group = "Physical data"));
        parameter Real shock_factor(final min = 0, start = 2.1) = 1
          "Shock factor dependent on cycle time and duty cycle (see diagram in Wittenstein catalogue)";
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = ratio) annotation(Placement(transformation(extent = {{0, -10}, {20, 10}})));
        BasicComponents.EfficiencyFactor efficiencyFactor(eta = eta, w_ref = 0.1 * n_1_N_lim * 2 * Constants.pi / 60) annotation(Placement(transformation(extent = {{-60, -10}, {-40, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Flange of right shaft"                                                          annotation(Placement(transformation(extent = {{78, -10}, {98, 10}}), iconTransformation(extent = {{78, -10}, {98, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        SI.Torque con_M_2 "Max acceleration torque at output";
        SI.Torque con_M_2_N
          "Maximum permitted output torque for continous duty";
        SI.Conversions.NonSIunits.AngularVelocity_rpm con_n_1
          "Maximum permitted input speed for short-time duty";
        SI.Conversions.NonSIunits.AngularVelocity_rpm con_n_1_N
          "Maximum permitted input speed for long-time duty";
        SI.Conversions.NonSIunits.AngularVelocity_rpm n_1_max
          "Maximum occurring velocity at input shaft";
        SI.Conversions.NonSIunits.AngularVelocity_rpm n_1_mean
          "Maximum occurring velocity at input shaft";
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J = J_GA) annotation(Placement(transformation(extent = {{-30, -10}, {-10, 10}})));
        BasicComponents.NoLoad_Friction noLoadFriction(M_Fr = M_noLoad_12) annotation(Placement(transformation(extent = {{-90, -10}, {-70, 10}})));
        Req_Wittenstein s annotation(Placement(transformation(extent = {{58, -10}, {78, 10}})));
        Modelica.Mechanics.Rotational.Components.Spring spring(c = c_t12) annotation(Placement(transformation(extent = {{30, -10}, {50, 10}})));
      equation
        n_1_max = s.n_max * ratio;
        n_1_mean = s.n_mean * ratio;
        con_M_2 = s.M_max * shock_factor - M_2_lim;
        con_M_2_N = s.M_eff3 - M_2_N_lim;
        con_n_1 = n_1_max - n_1_lim;
        con_n_1_N = n_1_mean - n_1_N_lim;
        connect(flange_a, noLoadFriction.flange_a) annotation(Line(points = {{-100, 0}, {-90, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(inertia.flange_a, efficiencyFactor.flange_b) annotation(Line(points = {{-30, 0}, {-40, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(noLoadFriction.flange_b, efficiencyFactor.flange_a) annotation(Line(points = {{-70, 0}, {-60, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(inertia.flange_b, idealGear.flange_a) annotation(Line(points = {{-10, 0}, {0, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(flange_b, s.flange_b) annotation(Line(points = {{88, 0}, {78, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGear.flange_b, spring.flange_a) annotation(Line(points = {{20, 0}, {30, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(s.flange_a, spring.flange_b) annotation(Line(points = {{58, 0}, {50, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-120, -100}, {120, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-120, -100}, {120, 100}}), graphics={  Rectangle(extent=  {{-100, 100}, {100, -100}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-100, 100}, {100, -100}}, fillColor=  {95, 95, 95},
                  fillPattern=                                                                                                    FillPattern.Solid, pattern=  LinePattern.None, lineColor=  {0, 0, 0}), Ellipse(extent=  {{-80, 80}, {80, -80}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-40, 40}, {40, -40}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-80, 20}, {-40, -20}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{40, 20}, {80, -20}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-53, -30}, {-13, -70}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{13, -30}, {53, -70}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-53, 70}, {-13, 30}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{13, 70}, {53, 30}}, lineColor=  {0, 0, 0}, fillColor=  {135, 135, 135},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-84, 24}, {-36, -24}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-56, 74}, {-8, 26}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{9, 73}, {57, 25}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{36, 24}, {84, -24}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{9, -26}, {57, -74}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-57, -26}, {-9, -74}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-36, 35}, {36, -36}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash), Ellipse(extent=  {{-84, 83}, {84, -84}}, lineColor=  {0, 0, 0}, pattern=  LinePattern.Dash)}));
      end Servogear_Wittenstein;

      model Req_Wittenstein
        Modelica.Mechanics.Rotational.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent = {{36, 70}, {56, 90}})));
        Modelica.Blocks.Interfaces.RealOutput M_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, 50}, {120, 70}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-96, 40}, {-76, 60}})));
        HelpBlocks.RootMeanSquareValue_custom rootMeanSquareValue_custom annotation(Placement(transformation(extent = {{58, -40}, {78, -20}})));
        Modelica.Blocks.Interfaces.RealOutput M_eff3
          "Connector of Real output signal"                                            annotation(Placement(transformation(extent = {{100, -40}, {120, -20}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Math.UnitConversions.To_rpm meanSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm absMaxSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, 20}, {0, 40}})));
        Modelica.Blocks.Interfaces.RealOutput n_max
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, 20}, {76, 40}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Interfaces.RealOutput n_mean
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, -10}, {76, 10}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        HelpBlocks.MeanValue meanSpeed
          annotation (Placement(transformation(extent={{-56,-10},{-36,10}})));
        HelpBlocks.AbsMax absMaxSpeed
          annotation (Placement(transformation(extent={{-56,20},{-36,40}})));
        HelpBlocks.AbsMax absMaxTorque
          annotation (Placement(transformation(extent={{58,50},{78,70}})));
      equation
        connect(torqueSensor.flange_b, flange_b) annotation(Line(points = {{56, 80}, {88, 80}, {88, 0}, {100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flange_a) annotation(Line(points = {{-96, 50}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.tau, rootMeanSquareValue_custom.u) annotation(Line(points = {{38, 69}, {38, -30}, {56, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, rootMeanSquareValue_custom.ref) annotation(Line(points = {{-75, 50}, {20, 50}, {20, -38}, {56, -38}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(meanSpeed_to_rpm.y, n_mean) annotation(Line(points = {{1, 0}, {66, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxSpeed_to_rpm.y, n_max) annotation(Line(points = {{1, 30}, {66, 30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom.y, M_eff3) annotation(Line(points = {{79, -30}, {110, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(torqueSensor.flange_a, flange_a) annotation(Line(points = {{36, 80}, {-96, 80}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(meanSpeed.y, meanSpeed_to_rpm.u) annotation (Line(
            points={{-35,0},{-22,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.y, absMaxSpeed_to_rpm.u) annotation (Line(
            points={{-35,30},{-22,30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-58,0},{-64,0},{-64,50},{20,50},{20,-38},{56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-58,30},{-64,30},{-64,50},{20,50},{20,-38},{56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.u, rootMeanSquareValue_custom.u) annotation (Line(
            points={{56,60},{38,60},{38,-30},{56,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.y, M_max) annotation (Line(
            points={{79,60},{110,60}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Text(extent=  {{-150, 73}, {150, 113}}, lineColor=  {0, 0, 255}, textString=  "Req. Wittenstein"), Ellipse(extent=  {{-70, 70}, {70, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 70}, {0, 40}}, color=  {0, 0, 0}), Line(points=  {{22.9, 32.8}, {40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-22.9, 32.8}, {-40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{37.6, 13.7}, {65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-37.6, 13.7}, {-65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{0, 0}, {9.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-0.48, 31.6}, {18, 26}, {18, 57.2}, {-0.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-5, 5}, {5, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, -100}, {-80, 0}}, color=  {0, 0, 127})}));
      end Req_Wittenstein;

      model Req_SEW
        Modelica.Mechanics.Rotational.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent = {{36, 70}, {56, 90}})));
        Modelica.Blocks.Interfaces.RealOutput M_a_max
          "Connector of Real output signal"                                             annotation(Placement(transformation(extent = {{100, 50}, {120, 70}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-96, 40}, {-76, 60}})));
        HelpBlocks.RootMeanSquareValue_custom rootMeanSquareValue_custom annotation(Placement(transformation(extent = {{58, -40}, {78, -20}})));
        HelpBlocks.RootMeanSquareValue_custom rootMeanSquareValue_custom1(root = 8) annotation(Placement(transformation(extent = {{58, -70}, {78, -50}})));
        HelpBlocks.RootMeanSquareValue_custom rootMeanSquareValue_custom2(root = 1.2) annotation(Placement(transformation(extent = {{58, -100}, {78, -80}})));
        Modelica.Blocks.Interfaces.RealOutput M_a_3
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, -40}, {120, -20}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Interfaces.RealOutput M_a_8
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, -70}, {120, -50}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Interfaces.RealOutput M_a_1_2
          "Connector of Real output signal"                                             annotation(Placement(transformation(extent = {{100, -100}, {120, -80}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Math.UnitConversions.To_rpm meanSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm absMaxSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, 20}, {0, 40}})));
        Modelica.Blocks.Interfaces.RealOutput n_a_max
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, 20}, {76, 40}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Interfaces.RealOutput n_a_mean
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, -10}, {76, 10}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        HelpBlocks.AbsMax absMaxTorque
          annotation (Placement(transformation(extent={{62,50},{82,70}})));
        HelpBlocks.AbsMax absMaxSpeed
          annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
        HelpBlocks.MeanValue meanSpeed
          annotation (Placement(transformation(extent={{-58,-10},{-38,10}})));
      equation
        connect(torqueSensor.flange_b, flange_b) annotation(Line(points = {{56, 80}, {88, 80}, {88, 0}, {100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flange_a) annotation(Line(points = {{-96, 50}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.tau, rootMeanSquareValue_custom.u) annotation(Line(points = {{38, 69}, {38, -30}, {56, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, rootMeanSquareValue_custom.ref) annotation(Line(points = {{-75, 50}, {20, 50}, {20, -38}, {56, -38}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom1.u, rootMeanSquareValue_custom.u) annotation(Line(points = {{56, -60}, {38, -60}, {38, -30}, {56, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom1.ref, rootMeanSquareValue_custom.ref) annotation(Line(points = {{56, -68}, {20, -68}, {20, -38}, {56, -38}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom2.u, rootMeanSquareValue_custom.u) annotation(Line(points = {{56, -90}, {38, -90}, {38, -30}, {56, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom2.ref, rootMeanSquareValue_custom.ref) annotation(Line(points = {{56, -98}, {20, -98}, {20, -38}, {56, -38}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(meanSpeed_to_rpm.y, n_a_mean) annotation(Line(points = {{1, 0}, {66, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxSpeed_to_rpm.y, n_a_max) annotation(Line(points = {{1, 30}, {66, 30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom.y, M_a_3) annotation(Line(points = {{79, -30}, {110, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(torqueSensor.flange_a, flange_a) annotation(Line(points = {{36, 80}, {-96, 80}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom2.y, M_a_1_2) annotation(Line(points = {{79, -90}, {110, -90}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom1.y, M_a_8) annotation(Line(points = {{79, -60}, {110, -60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxTorque.u, rootMeanSquareValue_custom.u) annotation (Line(
            points={{60,60},{38,60},{38,-30},{56,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.y, M_a_max) annotation (Line(
            points={{83,60},{110,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.y, absMaxSpeed_to_rpm.u) annotation (Line(
            points={{-37,30},{-22,30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-60,30},{-64,30},{-64,50},{20,50},{20,-38},{56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.y, meanSpeed_to_rpm.u) annotation (Line(
            points={{-37,0},{-22,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-60,0},{-64,0},{-64,50},{20,50},{20,-38},{56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Text(extent=  {{-150, 73}, {150, 113}}, lineColor=  {0, 0, 255}, textString=  "Req. SEW"), Ellipse(extent=  {{-70, 70}, {70, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 70}, {0, 40}}, color=  {0, 0, 0}), Line(points=  {{22.9, 32.8}, {40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-22.9, 32.8}, {-40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{37.6, 13.7}, {65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-37.6, 13.7}, {-65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{0, 0}, {9.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-0.48, 31.6}, {18, 26}, {18, 57.2}, {-0.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-5, 5}, {5, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, -100}, {-80, 0}}, color=  {0, 0, 127})}));
      end Req_SEW;
    end Servogear;

    package Clutch
      model Clutch_Wittenstein
        "Model of a clutch according to Wittenstein clutch catalogue"

        Modelica.Mechanics.Rotational.Components.Inertia J_1(J = J1) annotation(Placement(transformation(extent = {{-50, -10}, {-30, 10}})));
        Modelica.Mechanics.Rotational.Components.Inertia J_2(J = J2) annotation(Placement(transformation(extent = {{10, -10}, {30, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Req_Wittenstein s annotation(Placement(transformation(extent = {{40, -10}, {60, 10}})));
        Modelica.Mechanics.Rotational.Components.SpringDamper springDamper(c = c, d = d) annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        parameter SI.Inertia J1(min = 0, start = 0.001) "Moment of inertia" annotation(Dialog(group = "Physical data"));
        parameter SI.Inertia J2(min = 0, start = 0.001) "Moment of inertia" annotation(Dialog(group = "Physical data"));
        parameter SI.RotationalSpringConstant c(final min = 0, start = 1.0e5)
          "Dynamic torsional stiffness"                                                                     annotation(Dialog(group = "Physical data"));
        parameter Real psi = 1 "Proportional damping factor" annotation(Dialog(group = "Physical data"));
        parameter SI.Frequency f_test = 10 "Test Frequency, usually 10 Hz" annotation(Dialog(group = "Engineering data"));
        parameter SI.Torque M_lim = 310
          "Max acceleration torque (max 1000 cycles per hour)"                               annotation(Dialog(group = "Limiting data"));
        parameter SI.Torque M_N_lim = 130
          "Maximum permitted torque for short-time duty"                                 annotation(Dialog(group = "Limiting data"));
        parameter Real shock_factor(final min = 0, max = 2) = 1
          "Shock factor dependent on cycle time and duty cycle (see Wittenstein catlogue)";
        parameter Real temperature_factor(final min = 0, max = 2) = 1
          "Temperature factor (see Wittenstein catalogue)";
        SI.Torque con_M "Max acceleration torque";
        SI.Torque con_M_N "Maximum permitted torque for continous duty";
      protected
        parameter SI.RotationalDampingConstant d = psi * c / (2 * Constants.pi * omega)
          "Damping constant";
        parameter SI.AngularVelocity omega = 2 * Constants.pi * f_test;
      equation
        con_M = s.M_max * shock_factor * temperature_factor - M_lim;
        con_M_N = s.M_eff3 * shock_factor - M_N_lim;
        connect(flange_b, s.flange_b) annotation(Line(points = {{100, 0}, {60, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(flange_a, J_1.flange_a) annotation(Line(points = {{-100, 0}, {-50, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(springDamper.flange_a, J_1.flange_b) annotation(Line(points = {{-20, 0}, {-30, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(springDamper.flange_b, J_2.flange_a) annotation(Line(points = {{0, 0}, {10, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(s.flange_a, J_2.flange_b) annotation(Line(points = {{40, 0}, {30, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 40}, {-30, -40}}, lineColor=  {0, 0, 0}), Rectangle(extent=  {{30, 40}, {100, -40}}, lineColor=  {0, 0, 0}), Polygon(points=  {{-30, 40}, {-30, 60}, {-10, 60}, {-10, -20}, {6, -20}, {6, -60}, {-30, -60}, {-30, 40}}, lineColor=  {0, 0, 0}, smooth=  Smooth.None), Polygon(points=  {{-18, 40}, {-18, 60}, {2, 60}, {2, -20}, {18, -20}, {18, -60}, {-18, -60}, {-18, 40}}, lineColor=  {0, 0, 0}, smooth=  Smooth.None, origin=  {12, 0}, rotation=  180), Rectangle(extent=  {{-6, 20}, {6, -20}}, lineColor=  {0, 0, 0}), Text(extent=  {{-150, 59}, {150, 99}}, textString=  "%name", lineColor=  {0, 0, 255})}));
      end Clutch_Wittenstein;

      model Req_Wittenstein
        Modelica.Mechanics.Rotational.Sensors.TorqueSensor torqueSensor annotation(Placement(transformation(extent = {{36, 70}, {56, 90}})));
        Modelica.Blocks.Interfaces.RealOutput M_max
          "Connector of Real output signal"                                           annotation(Placement(transformation(extent = {{100, 50}, {120, 70}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
          "Left flange of shaft"                                                          annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
          "Right flange of shaft"                                                          annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(transformation(extent = {{-96, 40}, {-76, 60}})));
        HelpBlocks.RootMeanSquareValue_custom rootMeanSquareValue_custom annotation(Placement(transformation(extent = {{58, -40}, {78, -20}})));
        Modelica.Blocks.Interfaces.RealOutput M_eff3
          "Connector of Real output signal"                                            annotation(Placement(transformation(extent = {{100, -40}, {120, -20}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Math.UnitConversions.To_rpm meanSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, -10}, {0, 10}})));
        Modelica.Blocks.Math.UnitConversions.To_rpm absMaxSpeed_to_rpm annotation(Placement(transformation(extent = {{-20, 20}, {0, 40}})));
        Modelica.Blocks.Interfaces.RealOutput n_max
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, 20}, {76, 40}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        Modelica.Blocks.Interfaces.RealOutput n_mean
          "Connector of Real output signal containing input signal u in another unit"
                                                                                                              annotation(Placement(transformation(extent = {{56, -10}, {76, 10}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-80, -110})));
        HelpBlocks.AbsMax absMaxSpeed
          annotation (Placement(transformation(extent={{-54,20},{-34,40}})));
        HelpBlocks.MeanValue meanSpeed
          annotation (Placement(transformation(extent={{-54,-10},{-34,10}})));
        HelpBlocks.AbsMax absMaxTorque
          annotation (Placement(transformation(extent={{56,50},{76,70}})));
      equation
        connect(torqueSensor.flange_b, flange_b) annotation(Line(points = {{56, 80}, {88, 80}, {88, 0}, {100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(speedSensor.flange, flange_a) annotation(Line(points = {{-96, 50}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(torqueSensor.tau, rootMeanSquareValue_custom.u) annotation(Line(points = {{38, 69}, {38, -30}, {56, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(speedSensor.w, rootMeanSquareValue_custom.ref) annotation(Line(points = {{-75, 50}, {20, 50}, {20, -38}, {56, -38}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(meanSpeed_to_rpm.y, n_mean) annotation(Line(points = {{1, 0}, {66, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(absMaxSpeed_to_rpm.y, n_max) annotation(Line(points = {{1, 30}, {66, 30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(rootMeanSquareValue_custom.y, M_eff3) annotation(Line(points = {{79, -30}, {110, -30}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(torqueSensor.flange_a, flange_a) annotation(Line(points = {{36, 80}, {-96, 80}, {-96, 0}, {-100, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(absMaxSpeed.y, absMaxSpeed_to_rpm.u) annotation (Line(
            points={{-33,30},{-22,30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-56,30},{-62,30},{-62,50},{20,50},{20,-38},{56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.y, meanSpeed_to_rpm.u) annotation (Line(
            points={{-33,0},{-22,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(meanSpeed.u, rootMeanSquareValue_custom.ref) annotation (Line(
            points={{-56,0},{-60,0},{-60,30},{-62,30},{-62,50},{20,50},{20,-38},
                {56,-38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.y, M_max) annotation (Line(
            points={{77,60},{110,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(absMaxTorque.u, rootMeanSquareValue_custom.u) annotation (Line(
            points={{54,60},{38,60},{38,-30},{56,-30}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                  -100},{100,100}}),                                                                           graphics), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points=  {{-70, 0}, {-90, 0}}, color=  {0, 0, 0}), Line(points=  {{70, 0}, {90, 0}}, color=  {0, 0, 0}), Text(extent=  {{-150, 73}, {150, 113}}, lineColor=  {0, 0, 255}, textString=  "Req. Wittenstein"), Ellipse(extent=  {{-70, 70}, {70, -70}}, lineColor=  {0, 0, 0}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{0, 70}, {0, 40}}, color=  {0, 0, 0}), Line(points=  {{22.9, 32.8}, {40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{-22.9, 32.8}, {-40.2, 57.3}}, color=  {0, 0, 0}), Line(points=  {{37.6, 13.7}, {65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{-37.6, 13.7}, {-65.8, 23.9}}, color=  {0, 0, 0}), Line(points=  {{0, 0}, {9.02, 28.6}}, color=  {0, 0, 0}), Polygon(points=  {{-0.48, 31.6}, {18, 26}, {18, 57.2}, {-0.48, 31.6}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Ellipse(extent=  {{-5, 5}, {5, -5}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, -100}, {-80, 0}}, color=  {0, 0, 127})}));
      end Req_Wittenstein;
    end Clutch;
  end RotationalComponents;

  package Cooling
    package Interfaces
      connector HydraulicPort_b
        extends HydraulicPort;
        annotation(defaultComponentName = "port_b", Diagram(graphics={  Text(rotation = 0, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.None,
                  lineThickness =                                                                                                   0.25, extent = {{-150, 110}, {150, 50}}, textString = "%name"), Ellipse(rotation = 0, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}), Icon(graphics={  Ellipse(rotation = 0, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.Solid, extent = {{-40, 40}, {40, -40}})}));
      end HydraulicPort_b;

      connector HydraulicPort_a
        extends HydraulicPort;
        annotation(defaultComponentName = "port_a", Diagram(graphics={  Ellipse(rotation = 0, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.Solid,
                  lineThickness =                                                                                                   0.25, extent = {{-20, 20}, {20, -20}}), Text(rotation = 0, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.None,
                  lineThickness =                                                                                                   0.25, extent = {{-150, 110}, {150, 50}}, textString = "%name")}), Icon(graphics={  Ellipse(rotation = 0, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, pattern = LinePattern.Solid,
                  fillPattern =                                                                                                   FillPattern.Solid,
                  lineThickness =                                                                                                   0.25, extent = {{-40, 40}, {40, -40}})}));
      end HydraulicPort_a;

      connector HydraulicPort
        "Simplification the FluidPort from the Standard Library"
        flow SI.VolumeFlowRate q(displayUnit = "l/min")
          "Volume flow rate from the connection point into the component";
        Modelica.SIunits.Pressure p(start = 1, displayUnit = "bar")
          "Thermodynamic pressure in the connection point";
        Modelica.SIunits.Temperature T(start = 300, displayUnit = "degC");
      end HydraulicPort;

      model Water
        parameter SI.DynamicViscosity eta = 0.001 "Viscosity of the fluid";
        parameter SI.SpecificHeatCapacity c_p = 4182;
        parameter SI.Density rho = 1000;
      end Water;

      model TwoPort_Component
        replaceable model FluidProp = Water "constants of fluid models" annotation(Dialog(tab = "Fluid"));
        FluidProp fluidProp;
        HydraulicPort_a port_a annotation(Placement(transformation(extent = {{-114, -10}, {-94, 10}}), iconTransformation(extent = {{-114, -10}, {-94, 10}})));
        HydraulicPort_b port_b annotation(Placement(transformation(extent = {{94, -10}, {114, 10}}), iconTransformation(extent = {{94, -10}, {114, 10}})));
      end TwoPort_Component;
    end Interfaces;

    package Units
      type FlowResistanceT = Real(final quantity = "FlowResistance", final unit = "Pa.m-3.s");
      type FlowResistanceL = Real(final quantity = "FlowResistance", final unit = "Pa.m-3.s");
    end Units;

    model DisplacementPump
      "sets the flow rate by an input signal and defines temperature of the fluid"

      Interfaces.HydraulicPort_a port_a annotation(Placement(visible = true, transformation(origin = {0, -85.2829}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.SIunits.VolumeFlowRate Q(displayUnit = "l/min");
      parameter Modelica.SIunits.Temperature T "Temperature of liquid";
      Interfaces.HydraulicPort_b port_b annotation(Placement(transformation(extent = {{-10, 74}, {10, 94}}), iconTransformation(extent = {{-10, 70}, {10, 90}})));
      Modelica.Blocks.Interfaces.RealInput u annotation(Placement(transformation(extent = {{-92, -20}, {-52, 20}}), iconTransformation(extent = {{-92, -20}, {-52, 20}})));
    equation
      Q = u;
      port_b.q = -Q;
      port_b.q + port_a.q = 0;
      port_b.T = T;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent=  {{-60, 60}, {60, -60}}, lineColor=  {0, 0, 0}), Line(points=  {{0, 76}, {0, 60}}, color=  {0, 0, 0}, smooth=  Smooth.None), Line(points=  {{0, -80}, {0, -60}}, color=  {0, 0, 0}, smooth=  Smooth.None), Polygon(points=  {{20, 20}, {0, 60}, {-20, 20}, {20, 20}}, lineColor=  {0, 0, 0}, smooth=  Smooth.None, fillColor=  {0, 0, 0},
                fillPattern=                                                                                                    FillPattern.Solid)}));
    end DisplacementPump;

    model DisplacementPump_Constant
      "sets the flow rate by a parameter and defines temperature of the fluid"

      Interfaces.HydraulicPort_a port_a annotation(Placement(visible = true, transformation(origin = {0, -85.2829}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      parameter Modelica.SIunits.VolumeFlowRate Q(displayUnit = "l/min");
      parameter Modelica.SIunits.Temperature T "Temperature of liquid";
      Interfaces.HydraulicPort_b port_b annotation(Placement(transformation(extent = {{-10, 74}, {10, 94}}), iconTransformation(extent = {{-10, 70}, {10, 90}})));
    equation
      port_b.q = -Q;
      port_b.q + port_a.q = 0;
      port_b.T = T;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-60, 60}, {60, -60}}, lineColor = {0, 0, 0}), Line(points = {{0, 76}, {0, 60}}, color = {0, 0, 0}, smooth = Smooth.None), Line(points = {{0, -80}, {0, -60}}, color = {0, 0, 0}, smooth = Smooth.None), Polygon(points = {{20, 20}, {0, 60}, {-20, 20}, {20, 20}}, lineColor = {0, 0, 0}, smooth = Smooth.None, fillColor = {0, 0, 0},
                fillPattern =                                                                                                   FillPattern.Solid)}));
    end DisplacementPump_Constant;

    model Tank
      Interfaces.HydraulicPort_a port_a annotation(Placement(visible = true, transformation(origin = {0.15817, 60.348}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {-1.38778e-016, 60}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      parameter Modelica.SIunits.Pressure eps = 0.001
        "eps to avoid numerical problems";
    equation
      port_a.p = eps;
      annotation(Diagram, Icon(graphics={  Line(points = {{-40, 20}, {-40, -20}, {40.0351, -20.0351}, {40.0351, 20.3866}}, rotation = 0, color = {0, 0, 0}, pattern = LinePattern.Solid, thickness = 0.25), Line(points = {{0, 60}, {0, -20}}, rotation = 0, color = {0, 0, 0}, pattern = LinePattern.Solid, thickness = 0.25)}));
    end Tank;

    model CoolingChannels "Model for temperature and pressure drop"
      extends Interfaces.TwoPort_Component;
      parameter Units.FlowResistanceL R_flow "Resistance to flow of pipe";
      parameter SI.ThermalResistance R_therm "Thermal resistance";
      SI.Pressure delta_p(displayUnit = "bar") "Pressure drop at the channel";
      SI.MassFlowRate m_flow "Mass flow rate";
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport_a annotation(Placement(transformation(extent = {{-106, -10}, {-86, 10}}), iconTransformation(extent = {{-10, -120}, {10, -100}})));
      Modelica.SIunits.HeatFlowRate Q_flow_cooling;
    equation
      m_flow = port_a.q * fluidProp.rho;
      port_b.T = heatport_a.T + (port_a.T - heatport_a.T) * exp(-1 / (R_therm * m_flow * fluidProp.c_p));
      Q_flow_cooling = port_a.q * fluidProp.rho * fluidProp.c_p * (port_b.T - port_a.T);
      Q_flow_cooling = heatport_a.Q_flow;
      delta_p = sign(port_a.q) * R_flow * abs(port_a.q) ^ 1.75; //Based on the law of Blasius
      port_a.q + port_b.q = 0;
      delta_p = port_a.p - port_b.p;
      annotation(Documentation(info = "<html>
<p>
Assumptions:
<ul>
<li>Fluid is viscous and incompressible</li>
<li>Flow is laminar through a pipe of constant circular cross-section</li>
<li>Length >> Diameter</li>
<li>No acceleration of fluid in the pipe</li>
</ul>
</p>
</HTML>
"), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent=  {{-100, 60}, {100, 40}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                fillPattern=                                                                                                    FillPattern.Forward), Rectangle(extent=  {{-100, -40}, {100, -60}}, lineColor=  {0, 0, 0}, fillColor=  {0, 0, 0},
                fillPattern=                                                                                                    FillPattern.Forward), Line(points=  {{-100, 0}, {-74, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{-70, 0}, {-68, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{-64, 0}, {-38, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{-34, 0}, {-32, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{-28, 0}, {-2, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{2, 0}, {4, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{8, 0}, {34, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{38, 0}, {40, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{44, 0}, {70, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{74, 0}, {76, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Line(points=  {{80, 0}, {100, 0}}, color=  {0, 0, 0}, thickness=  0.5, smooth=  Smooth.None), Text(extent=  {{100, 0}, {-100, 40}}, lineColor=  {0, 0, 0},
                lineThickness=                                                                                                    0.5, fillColor=  {0, 0, 0},
                fillPattern=                                                                                                    FillPattern.Forward, textString=  "flow"), Polygon(points=  {{-10, -4}, {-10, 15}, {0, 35}, {10, 15}, {10, -4}, {-10, -4}}, lineColor=  {255, 0, 0}, origin=  {0, -96}, rotation=  360), Text(extent=  {{-150, 125}, {150, 85}}, textString=  "%name", lineColor=  {0, 0, 255})}));
    end CoolingChannels;
  end Cooling;

  package HelpBlocks

    block Max
      "Computes the maximum of a signal. Similar to max block in the optimization library but without differentiation"
      extends Modelica.Blocks.Interfaces.SISO;
      Real y_delay(start = 0);
    protected
      constant Real eps = 1e-4;
    equation
      y_delay = delay(y, eps);
    algorithm
      if time < eps then
        y := u;
      else
        y := max(u, y_delay);
      end if;
      annotation(Icon(coordinateSystem(preserveAspectRatio=false,  extent={{-100,
                -100},{100,100}},                                                                       grid = {1, 1}), graphics={Text(
              extent={{100,-20},{-100,20}},
              lineColor={0,0,0},
              textString="max(f(t))")}),                                                                                                    Documentation(info = "<html>
<p>The <b>maximum</b> of the input signal u is detected.</p>
</html>"));
    end Max;

    block MeanValue "Mean value of a signal"
      extends Modelica.Blocks.Interfaces.SISO;
    initial equation
      y = 0;
    equation
      der(y*time) = u;
      annotation(Icon(coordinateSystem(preserveAspectRatio=false,  extent={{-100,-100},
                {100,100}},                                                                             grid = {1, 1}), graphics={Text(
              extent={{-100,20},{100,-20}},
              lineColor={0,0,0},
              textString="mean")}));
    end MeanValue;

    block AbsMax "Maximum absolute value of a input"
      extends Modelica.Blocks.Interfaces.SISO;
      Modelica.Blocks.Math.Abs abs1 annotation(Placement(transformation(extent = {{-40, -10}, {-20, 10}})));
      Max max annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    equation
      connect(u, abs1.u) annotation(Line(points = {{-120, 0}, {-42, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(abs1.y, max.u) annotation (Line(
          points={{-19,0},{18,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(max.y, y) annotation (Line(
          points={{41,0},{110,0}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={                                                                                                    Rectangle(extent = {{-80, -9}, {-5, -41}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid),                                                                                                    Rectangle(extent = {{-80, -9}, {-54, -24}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-12, -9}, {82, -63}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid),
                                                                                                    Text(
              extent={{100,-20},{-100,20}},
              lineColor={0,0,0},
              textString="Abs
[max(f(t))]")}),
           Diagram(coordinateSystem(preserveAspectRatio=false,  extent={{-100,-100},
                {100,100}}),                                                                          graphics));
    end AbsMax;

    block DynOpt "Computes the objectives for control optimization"

      Real obj_meanDeviation;
      Real obj_maxOvershoot;
      Real obj_settlingTime;
      Real obj_risingTime;
      Real con_maxOvershoot;
      Real allowedOvershoot;
      Real obj_ITAE( start=0);
      Real temp;
      parameter Real stepHeight = 1 "hight of input step";
      parameter Real allowedOvershoot_percent = 0 "% of step height";
      parameter Real settlingTolerance_percent = 3 "% of step height";

      Modelica.Blocks.Interfaces.RealInput ref annotation (Placement(transformation(
              extent={{-140,-100},{-100,-60}}),
                                             iconTransformation(extent={{-140,-100},
                {-100,-60}})));
      Modelica.Blocks.Interfaces.RealInput act annotation (Placement(transformation(
              extent={{-140,60},{-100,100}}),  iconTransformation(extent={{-140,60},
                {-100,100}})));
      Feedback_mirror            pos_feedback annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-48,-34})));
      Modelica.Blocks.Math.Abs abs1
        annotation (Placement(transformation(extent={{2,30},{22,50}})));
      HelpBlocks.Max max
        annotation (Placement(transformation(extent={{2,70},{22,90}})));
      HelpBlocks.MeanValue meanValue
        annotation (Placement(transformation(extent={{42,30},{62,50}})));
      HelpBlocks.Max max_actMinusRef
        annotation (Placement(transformation(extent={{8,-44},{28,-24}})));
    initial equation
      obj_settlingTime = 0;
      temp = 0;
    equation
      allowedOvershoot = stepHeight*( 1+ allowedOvershoot_percent / 100);
      obj_meanDeviation = meanValue.y;
      obj_maxOvershoot = max.y;
      der(obj_ITAE) = time*abs(act-ref);
      con_maxOvershoot = max.y - allowedOvershoot;
      when {abs(ref-act) < settlingTolerance_percent/100*stepHeight} then
        obj_settlingTime = time;
      end when;

      connect(pos_feedback.u1, act) annotation (Line(
          points={{-56,-34},{-66,-34},{-66,80},{-120,80}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ref, pos_feedback.u2) annotation (Line(
          points={{-120,-80},{-48,-80},{-48,-26}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pos_feedback.y, abs1.u) annotation (Line(
          points={{-39,-34},{-26,-34},{-26,40},{0,40}},
          color={0,0,127},
          smooth=Smooth.None));

    algorithm
      when {abs(ref-act) < settlingTolerance_percent/100*stepHeight and temp == 0} then
        obj_risingTime :=time;
        temp :=1;
      end when;

    equation
      connect(max.u, act) annotation (Line(
          points={{0,80},{-120,80}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(abs1.y, meanValue.u) annotation (Line(
          points={{23,40},{40,40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(max_actMinusRef.u, abs1.u) annotation (Line(
          points={{6,-34},{-26,-34},{-26,40},{0,40}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Icon(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                -100},{100,100}},                                                                        grid = {1, 1}), graphics={
                                    Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
            Text(
              extent={{-192,-100},{-100,-150}},
              lineColor={0,0,255},
              textString="ref"),                                                                                                    Polygon(points={{
                  -56,67},{-64,45},{-48,45},{-56,67}},                                                                                                    lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points={{
                  -56,59},{-56,-61}},                                                                                                    color = {192, 192, 192}), Line(points={{
                  -56,-61},{43,-61}},                                                                                                    color = {192, 192, 192}), Polygon(points={{
                  64,-60},{42,-52},{42,-68},{64,-60}},                                                                                                    lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid),
            Line(
              points={{-55,-61},{-18,-61}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-18,-61},{-18,5},{-18,39}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-18,39},{58,39}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-18,-62},{-7,-27},{12,76},{36,29},{54,47}},
              color={0,0,0},
              smooth=Smooth.Bezier),        Text(
            extent={{-150,160},{150,120}},
            textString="%name",
            lineColor={0,0,255}),
            Text(
              extent={{-192,151},{-100,101}},
              lineColor={0,0,255},
              textString="act")}),                                                                                                    Documentation(info = "<html>
<p>The<b> mean value</b> y is an integral time-averaged value of the input variable u. </p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                        graphics));
    end DynOpt;

    block RootMeanSquareValue
      "Computes the integral time-averaged root mean value of a signal"
      extends Modelica.Blocks.Interfaces.SISO;
    protected
      Real y1;
      parameter Real startTime(fixed = false);
    initial equation
      startTime = time;
      y1 = 0;
      y = u;
    equation
      der(y1) = u ^ 2;
      when terminal() then
        y = sqrt(y1 / (time - startTime));
      end when;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Text(extent = {{-100, 41}, {100, -40}}, lineColor = {0, 0, 255}, fillColor = {255, 170, 85},
                fillPattern =                                                                                                   FillPattern.Solid, textString = "RMS")}), Documentation(info = "<html>
<p>The<b> mean value</b> y is an integral time-averaged value of the input variable u. </p>
</html>"));
    end RootMeanSquareValue;

    block RootMeanSquareValue_custom
      "Computes the weighted integral time-averaged mean value of a signal by using the p-norm"
      extends Modelica.Blocks.Interfaces.SISO;
      parameter Real root = 3 "Choose e.g. 3 for cubic root";
    protected
      Real numerator;
      Real denominator;
      parameter Real startTime(fixed = false);
    public
      Modelica.Blocks.Interfaces.RealInput ref annotation(Placement(transformation(extent = {{-140, -100}, {-100, -60}}), iconTransformation(extent = {{-140, -100}, {-100, -60}})));
    initial equation
      startTime = time;
      numerator = 0;
      denominator = 0;
      y = u;
    equation
      der(numerator) = abs(ref) * abs(u) ^ root;
      der(denominator) = abs(ref);
      when terminal() then
        y = (numerator / denominator) ^ (1 / root);
      end when;
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Text(extent = {{-100, 41}, {100, -40}}, lineColor = {0, 0, 255}, fillColor = {255, 170, 85},
                fillPattern =                                                                                                   FillPattern.Solid, textString = "RMS"), Text(extent = {{-100, 100}, {100, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.Dash, fillColor = {0, 0, 255},
                fillPattern =                                                                                                   FillPattern.Solid, textString = "root = %root"), Text(extent = {{-73, -55}, {-3, -101}}, lineColor = {0, 0, 255}, pattern = LinePattern.Dash, fillColor = {0, 0, 255},
                fillPattern =                                                                                                   FillPattern.Solid, textString = "ref",
                horizontalAlignment =                                                                                                   TextAlignment.Left)}), Documentation(info = "<html>
<p>The<b> mean value</b> y is an integral time-averaged value of the input variable u. </p>
</html>"));
    end RootMeanSquareValue_custom;

    block TimeTable_Periodic
      "Generate a (possibly discontinuous) signal by linear interpolation in a table and repeats the signal each periode"
      parameter Real table[:, 2]
        "Table matrix (time = first column; e.g., table=[0, 0; 1, 1; 2, 4])";
      parameter Real offset = 0 "Offset of output signal";
      parameter SI.Time startTime = 0 "Output = offset for time < startTime";
      extends Modelica.Blocks.Interfaces.SO;
      parameter Modelica.SIunits.Time period(final min = Modelica.Constants.small, start = 1)
        "Time for one period";
    protected
      Real a "Interpolation coefficients a of actual interval (y=a*x+b)";
      Real b "Interpolation coefficients b of actual interval (y=a*x+b)";
      Integer last(start = 1) "Last used lower grid index";
      SI.Time nextEvent(start = 0, fixed = true) "Next event instant";
      Modelica.SIunits.Time T_start(start = 0) "Start time of current period";
      Integer count(start=0, fixed = true) "Period count";

      function getInterpolationCoefficients
        "Determine interpolation coefficients and next time event"
        input Real table[:, 2] "Table for interpolation";
        input Real offset "y-offset";
        input Real startTime "time-offset";
        input Real t "Actual time instant";
        input Integer last "Last used lower grid index";
        input Real TimeEps
          "Relative epsilon to check for identical time instants";
        output Real a "Interpolation coefficients a (y=a*x + b)";
        output Real b "Interpolation coefficients b (y=a*x + b)";
        output Real nextEvent "Next event instant";
        output Integer next "New lower grid index";
      protected
        Integer columns = 2 "Column to be interpolated";
        Integer ncol = 2 "Number of columns to be interpolated";
        Integer nrow = size(table, 1) "Number of table rows";
        Integer next0;
        Real tp;
        Real dt;
      algorithm
        next := last;
        nextEvent := t - TimeEps * abs(t);
        // in case there are no more time events
        tp := t + TimeEps * abs(t) - startTime;
        if tp < 0.0 then
          nextEvent := startTime;
          a := 0;
          b := offset;
        elseif nrow < 2 then
          a := 0;
          b := offset + table[1, columns];
        else
          while next < nrow and tp >= table[next, 1] loop
            next := next + 1;
          end while;
          if next < nrow then
            nextEvent := startTime + table[next, 1];
          end if;
          next0 := next - 1;
          dt := table[next, 1] - table[next0, 1];
          if dt <= TimeEps * abs(table[next, 1]) then
            a := 0;
            b := offset + table[next, columns];
          else
            a := (table[next, columns] - table[next0, columns]) / dt;
            b := offset + table[next0, columns] - a * table[next0, 1];
          end if;
        end if;
        // First event not yet reached
        // Special action if table has only one row
        // Find next time event instant. Note, that two consecutive time instants
        // in the table may be identical due to a discontinuous point.
        // Define next time event, if last table entry not reached
        // Determine interpolation coefficients
        // Interpolation interval is not big enough, use "next" value
        // Take into account startTime "a*(time - startTime) + b"
        b := b - a * startTime;
      end getInterpolationCoefficients;
    algorithm
      when {time - T_start >= pre(nextEvent), initial()} then
        (a, b, nextEvent, last) := getInterpolationCoefficients(table, offset, startTime, time - T_start, last, 100 * Modelica.Constants.eps);
      end when;
      when integer((time - startTime) / period) > pre(count) then
        count := pre(count) + 1;
        T_start := time;
        nextEvent := 0;
        last := 1;
        (a, b, nextEvent, last) := getInterpolationCoefficients(table, offset, startTime, 0, last, 100 * Modelica.Constants.eps);
      end when;
    equation
      y = a * (time - T_start) + b;
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-48, 70}, {2, -50}}, lineColor = {255, 255, 255}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-48, -50}, {-48, 70}, {52, 70}, {52, -50}, {-48, -50}, {-48, -20}, {52, -20}, {52, 10}, {-48, 10}, {-48, 40}, {52, 40}, {52, 70}, {2, 70}, {2, -51}}, color = {0, 0, 0}), Text(extent = {{-150, -150}, {150, -110}}, lineColor = {0, 0, 0}, textString = "offset=%offset")}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics={  Polygon(points=  {{-80, 90}, {-85, 68}, {-74, 68}, {-80, 90}}, lineColor=  {95, 95, 95}, fillColor=  {95, 95, 95},
                fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, 68}, {-80, -80}}, color=  {95, 95, 95}), Line(points=  {{-90, -70}, {82, -70}}, color=  {95, 95, 95}), Polygon(points=  {{88, -70}, {68, -65}, {68, -74}, {88, -70}}, lineColor=  {95, 95, 95}, fillColor=  {95, 95, 95},
                fillPattern=                                                                                                    FillPattern.Solid), Rectangle(extent=  {{-20, 90}, {30, -30}}, lineColor=  {255, 255, 255}, fillColor=  {192, 192, 192},
                fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-20, -30}, {-20, 90}, {80, 90}, {80, -30}, {-20, -30}, {-20, 0}, {80, 0}, {80, 30}, {-20, 30}, {-20, 60}, {80, 60}, {80, 90}, {30, 90}, {30, -31}}, color=  {0, 0, 0}), Text(extent=  {{-70, -42}, {-32, -54}}, lineColor=  {0, 0, 0}, textString=  "offset"), Polygon(points=  {{-31, -30}, {-33, -40}, {-28, -40}, {-31, -30}}, lineColor=  {95, 95, 95}, fillColor=  {95, 95, 95},
                fillPattern=                                                                                                    FillPattern.Solid), Polygon(points=  {{-31, -70}, {-34, -60}, {-29, -60}, {-31, -70}, {-31, -70}}, lineColor=  {95, 95, 95}, fillColor=  {95, 95, 95},
                fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-31, -32}, {-31, -70}}, color=  {95, 95, 95}), Line(points=  {{-20, -30}, {-20, -70}}, color=  {95, 95, 95}), Text(extent=  {{-38, -73}, {8, -83}}, lineColor=  {0, 0, 0}, textString=  "startTime"), Line(points=  {{-20, -30}, {-80, -30}}, color=  {95, 95, 95}), Text(extent=  {{-76, 93}, {-44, 75}}, lineColor=  {0, 0, 0}, textString=  "y"), Text(extent=  {{66, -78}, {90, -88}}, lineColor=  {0, 0, 0}, textString=  "time"), Text(extent=  {{-15, 83}, {24, 68}}, lineColor=  {0, 0, 0}, textString=  "time"), Text(extent=  {{33, 83}, {76, 67}}, lineColor=  {0, 0, 0}, textString=  "y")}), Documentation(info = "<HTML>
<p>
This block generates an output signal by <b>linear interpolation</b> in
a table. The time points and function values are stored in a matrix
<b>table[i,j]</b>, where the first column table[:,1] contains the
time points and the second column contains the data to be interpolated.
The table interpolation has the following proporties:
</p>
<ul>
<li>The time points need to be <b>monotonically increasing</b>. </li>
<li><b>Discontinuities</b> are allowed, by providing the same
    time point twice in the table. </li>
<li>Values <b>outside</b> of the table range, are computed by
    <b>extrapolation</b> through the last or first two points of the
    table.</li>
<li>If the table has only <b>one row</b>, no interpolation is performed and
    the function value is just returned independantly of the
    actual time instant.</li>
<li>Via parameters <b>startTime</b> and <b>offset</b> the curve defined
    by the table can be shifted both in time and in the ordinate value.
<li>The table is implemented in a numerically sound way by
    generating <b>time events</b> at interval boundaries,
    in order to not integrate over a discontinuous or not differentiable
    points.
</li>
</ul>
<p>
Example:
</p>
<pre>
   table = [0  0
            1  0
            1  1
            2  4
            3  9
            4 16]
If, e.g., time = 1.0, the output y =  0.0 (before event), 1.0 (after event)
    e.g., time = 1.5, the output y =  2.5,
    e.g., time = 2.0, the output y =  4.0,
    e.g., time = 5.0, the output y = 23.0 (i.e., extrapolation).
</pre>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/TimeTable.png\">
</p>

</HTML>
",     revisions=
               "<html>
<p><b>Release Notes:</b></p>
<ul>
<li><i>Oct. 21, 2002</i>
       by <a href=\"http://www.robotic.dlr.de/Christian.Schweiger/\">Christian Schweiger</a>:<br>
       Corrected interface from
<pre>
    parameter Real table[:, :]=[0, 0; 1, 1; 2, 4];
</pre>
       to
<pre>
    parameter Real table[:, <b>2</b>]=[0, 0; 1, 1; 2, 4];
</pre>
       </li>
<li><i>Nov. 7, 1999</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"));
    end TimeTable_Periodic;

    model PolynomialFunction
      "Output the polynomial function of the absolute input"
      extends Modelica.Blocks.Interfaces.SISO;
      parameter Real exponent;
      parameter Real factor;
    equation
      y = factor * abs(u) ^ exponent;
      annotation(defaultComponentName = "sqrt1", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Line(points = {{-90, -80}, {68, -80}}, color = {192, 192, 192}), Polygon(points = {{90, -80}, {68, -72}, {68, -88}, {90, -80}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-80, -80}, {-79.2, -68.7}, {-78.4, -64}, {-76.8, -57.3}, {-73.6, -47.9}, {-67.9, -36.1}, {-59.1, -22.2}, {-46.2, -6.49}, {-28.5, 10.7}, {-4.42, 30}, {27.7, 51.3}, {69.5, 74.7}, {80, 80}}, color = {0, 0, 0}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-80, -88}, {-80, 68}}, color = {192, 192, 192}), Text(extent = {{-8, -4}, {64, -52}}, lineColor = {192, 192, 192}, textString = "poly")}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Line(points=  {{-92, -80}, {84, -80}}, color=  {192, 192, 192}), Polygon(points=  {{100, -80}, {84, -74}, {84, -86}, {100, -80}}, lineColor=  {192, 192, 192}, fillColor=  {192, 192, 192},
                fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, -80}, {-79.2, -68.7}, {-78.4, -64}, {-76.8, -57.3}, {-73.6, -47.9}, {-67.9, -36.1}, {-59.1, -22.2}, {-46.2, -6.49}, {-28.5, 10.7}, {-4.42, 30}, {27.7, 51.3}, {69.5, 74.7}, {80, 80}}, color=  {0, 0, 0}), Polygon(points=  {{-80, 98}, {-86, 82}, {-74, 82}, {-80, 98}}, lineColor=  {192, 192, 192}, fillColor=  {192, 192, 192},
                fillPattern=                                                                                                    FillPattern.Solid), Line(points=  {{-80, -90}, {-80, 84}}, color=  {192, 192, 192}), Text(extent=  {{-71, 98}, {-44, 78}}, lineColor=  {160, 160, 164}, textString=  "y"), Text(extent=  {{60, -52}, {84, -72}}, lineColor=  {160, 160, 164}, textString=  "u")}), Documentation(info = "<HTML>


</HTML>
"));
    end PolynomialFunction;

    model TemperatureSensor "computes acutual and maximum temperature"

      Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation(Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
      Modelica.Blocks.Interfaces.RealOutput max annotation(Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0), iconTransformation(extent = {{96, -20}, {116, 0}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0), iconTransformation(extent = {{-104, -10}, {-84, 10}})));
      Modelica.Blocks.Interfaces.RealOutput actual annotation(Placement(transformation(extent = {{92, -40}, {112, -20}}, rotation = 0), iconTransformation(extent = {{96, 0}, {116, 20}})));
      HelpBlocks.AbsMax absMaxTemperature
        annotation (Placement(transformation(extent={{18,-10},{38,10}})));
    equation
      connect(temperatureSensor.port, port) annotation(Line(points = {{-70, 0}, {-100, 0}}, color = {191, 0, 0}, smooth = Smooth.None));
      connect(temperatureSensor.T, absMaxTemperature.u) annotation (Line(
          points={{-50,0},{16,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(absMaxTemperature.y, max) annotation (Line(
          points={{39,0},{100,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(actual, absMaxTemperature.u) annotation (Line(
          points={{102,-30},{2,-30},{2,0},{16,0}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-14, -98}, {26, -60}}, lineColor = {0, 0, 0},
                lineThickness =                                                                                                   0.5, fillColor = {191, 0, 0},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-6, 40}, {18, -68}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{18, 0}, {96, 0}}, color = {0, 0, 255}), Line(points = {{-84, 0}, {-6, 0}}, color = {191, 0, 0}), Polygon(points = {{-6, 40}, {-6, 80}, {-4, 86}, {0, 88}, {6, 90}, {12, 88}, {16, 86}, {18, 80}, {18, 40}, {-6, 40}}, lineColor = {0, 0, 0},
                lineThickness =                                                                                                   0.5), Line(points = {{-6, 40}, {-6, -64}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{18, 40}, {18, -64}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-34, -20}, {-6, -20}}, color = {0, 0, 0}), Line(points = {{-34, 20}, {-6, 20}}, color = {0, 0, 0}), Line(points = {{-34, 60}, {-6, 60}}, color = {0, 0, 0}), Text(extent = {{132, -20}, {32, -120}}, lineColor = {0, 0, 0}, textString = "?C"), Text(extent = {{-144, 130}, {156, 90}}, textString = "%name", lineColor = {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio=false,   extent={{-100,
                -100},{100,100}}),                                                                                                    graphics));
    end TemperatureSensor;

    block Feedback_mirror
      "Mirror image of the standard feedbackOutput difference between commanded and feedback input"

      input Modelica.Blocks.Interfaces.RealInput u1 annotation(Placement(transformation(extent = {{-100, -20}, {-60, 20}}, rotation = 0)));
      input Modelica.Blocks.Interfaces.RealInput u2 annotation(Placement(transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
      output Modelica.Blocks.Interfaces.RealOutput y annotation(Placement(transformation(extent = {{80, -10}, {100, 10}}, rotation = 0)));
    equation
      y = u1 - u2;
      annotation(Documentation(info = "
<HTML>
<p>
This blocks computes output <b>y</b> as <i>difference</i> of the
commanded input <b>u1</b> and the feedback
input <b>u2</b>:
</p>
<pre>
    <b>y</b> = <b>u1</b> - <b>u2</b>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   n = 2

  results in the following equations:

     y = u1 - u2
</pre>

</HTML>
"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Ellipse(extent = {{-20, 20}, {20, -20}}, lineColor = {0, 0, 127}, fillColor = {235, 235, 235},
                fillPattern =                                                                                                   FillPattern.Solid), Line(points = {{-60, 0}, {-20, 0}}, color = {0, 0, 127}), Line(points = {{20, 0}, {80, 0}}, color = {0, 0, 127}), Line(points = {{0, -20}, {0, -60}}, color = {0, 0, 127}), Text(extent = {{-82, 2}, {14, -92}}, lineColor = {0, 0, 0}, textString = "-"), Text(extent = {{-150, 94}, {150, 44}}, textString = "%name", lineColor = {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Ellipse(extent=  {{-20, 20}, {20, -20}}, pattern=  LinePattern.Solid,
                lineThickness=                                                                                                    0.25, fillColor=  {235, 235, 235},
                fillPattern=                                                                                                    FillPattern.Solid, lineColor=  {0, 0, 255}), Line(points=  {{-60, 0}, {-20, 0}}, color=  {0, 0, 255}), Line(points=  {{20, 0}, {80, 0}}, color=  {0, 0, 255}), Line(points=  {{0, -20}, {0, -60}}, color=  {0, 0, 255}), Text(extent=  {{-86, 10}, {10, -84}}, lineColor=  {0, 0, 0}, textString=  "-")}));
    end Feedback_mirror;

    model Machining
      "defines axial and normal force caused by machining as well as the velocity profile"

      TranslationalComponents.Interfaces.Flange_b flange_b annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      TranslationalComponents.Interfaces.FlangeChangeLeft flangeChangeLeft annotation(Placement(transformation(extent = {{-68, -10}, {-48, 10}})));
      Modelica.Mechanics.Translational.Sources.Force force annotation(Placement(transformation(extent = {{60, 70}, {80, 90}})));
      Modelica.Mechanics.Translational.Sources.Force force1 annotation(Placement(transformation(extent = {{40, -50}, {60, -30}})));
      Modelica.Blocks.Sources.TimeTable tableForce_ax(table = force_ax) annotation(Placement(transformation(extent = {{-10, 70}, {10, 90}})));
      Modelica.Blocks.Sources.TimeTable tableForce_N(table = force_N) annotation(Placement(transformation(extent = {{-10, -50}, {10, -30}})));
      parameter Real force_ax[:, 2] = [0, 0] "Time table for axial force";
      parameter Real vel_ax[:, 2] = [0, 0] "Time table for axial velocity";
      parameter Real force_N[:, 2] = [0, 0] "Time table for normal force";
    equation
      connect(flangeChangeLeft.flange_a, flange_b) annotation(Line(points = {{-63, 0}, {-84, 0}, {-84, 0}, {-100, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(force.flange, flangeChangeLeft.flange_ax) annotation(Line(points = {{80, 80}, {80, 10}, {-53, 10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(force1.flange, flangeChangeLeft.flange_N) annotation(Line(points = {{60, -40}, {80, -40}, {80, -10}, {-53, -10}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(force.f, tableForce_ax.y) annotation(Line(points = {{58, 80}, {11, 80}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(tableForce_N.y, force1.f) annotation(Line(points = {{11, -40}, {38, -40}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Icon(graphics={  Rectangle(extent = {{-100, 62}, {100, -60}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-60, 0}, {60, -40}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{-2, 60}, {70, -12}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{10, 48}, {58, 0}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid, pattern = LinePattern.Dash), Text(extent = {{-150, 100}, {150, 60}}, textString = "%name", lineColor = {0, 0, 255})}));
    end Machining;




    function print
      input String str;
      input String file;
    algorithm
      Modelica.Utilities.Streams.print(str, file);
    end print;

    function printRealArray "Print string to terminal or file"
      extends Modelica.Icons.Function;
      input Real[:] x "Input to be printed";
      input String fileName = ""
        "File where to print (empty string is the terminal)";
      input Integer minimumLength = 1 "Minimum width of result";
      input Integer significantDigits = 6 "Number of significant digits";
      output String outStr = "" "String to be printed";
    algorithm
      for i in 1:size(x, 1) loop
        outStr := outStr + " " + String(x[i], minimumLength=  minimumLength, significantDigits=  significantDigits);
      end for;
      Modelica.Utilities.Streams.print(string=  outStr, fileName=  fileName);
    end printRealArray;

    block To_mmPerMin "Convert from metre per second to mm per minute"
      extends Modelica.Blocks.Interfaces.PartialConversionBlock(u(unit = "m/s"), y(unit = "mm/min"));
    equation
      y = 1000 * 60 * u;
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{0, 82}, {-96, 42}}, lineColor = {0, 0, 0}, textString = "m/s"), Text(extent = {{92, -40}, {-14, -84}}, lineColor = {0, 0, 0}, textString = "mm/min")}), Documentation(info = "<html>
<p>
This block converts the input signal from metre per second to kilometre per hour and returns
the result as output signal.
</p>
</html>"));
    end To_mmPerMin;

    block PositionController
      "Controller with unit conversion. Uses the units that are typical for feed drives."
      parameter Real ratio_R2T(start = 1, unit = "m") "mm per revolution";
      parameter Real K_v(start = 1, unit = "m/(mm.min)")
        "Gain value multiplied with input signal";
    public
      Modelica.Blocks.Interfaces.RealInput u "Input signal connector" annotation(Placement(transformation(extent = {{-100, -20}, {-60, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput y "Output signal connector" annotation(Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
    equation
      y = K_v / 60 / ratio_R2T * (2 * Constants.pi) * u;
      annotation(Documentation(info = "<html>
<p>
This block computes output <i>y</i> as
<i>product</i> of gain <i>k</i> with the
input <i>u</i>:
</p>
<pre>
    y = k * u;
</pre>

</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Polygon(points = {{-60, -60}, {-60, 60}, {60, 0}, {-60, -60}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                   FillPattern.Solid), Text(extent = {{-150, -140}, {150, -100}}, lineColor = {0, 0, 0}, textString = "K_v=%K_v"), Text(extent = {{-150, 140}, {150, 100}}, textString = "%name", lineColor = {0, 0, 255})}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Polygon(points=  {{-60, -60}, {-60, 60}, {60, 0}, {-60, -60}}, lineColor=  {0, 0, 127}, fillColor=  {255, 255, 255},
                fillPattern=                                                                                                    FillPattern.Solid), Text(extent=  {{-76, 38}, {0, -34}}, textString=  "k", lineColor=  {0, 0, 255})}));
    end PositionController;

    block from_radPerSec_to_mPerMin
      extends Modelica.Blocks.Interfaces.PartialConversionBlock(u(unit = "rad/s"), y(unit = "m/min"));
      parameter Real ratio_R2T(start = 1, unit = "m") "mm per revolution";
    equation
      y = 60 / (2 * Constants.pi) * ratio_R2T * u;
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{0, 82}, {-96, 42}}, lineColor = {0, 0, 0}, textString = "rad/s"), Text(extent = {{92, -40}, {-14, -84}}, lineColor = {0, 0, 0}, textString = "m/min")}), Documentation(info = "<html>
<p>
This block converts the input signal from metre per second to kilometre per hour and returns
the result as output signal.
</p>
</html>"));
    end from_radPerSec_to_mPerMin;

    block from_radPerSec_to_mmPerMin
      extends Modelica.Blocks.Interfaces.PartialConversionBlock(u(unit = "rad/s"), y(unit = "m/min"));
      parameter Real ratio_R2T(start = 1, unit = "mm") "mm per revolution";
    equation
      y = 60 / (2 * Constants.pi) * ratio_R2T * u;
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{0, 82}, {-96, 42}}, lineColor = {0, 0, 0}, textString = "rad/s"), Text(extent = {{92, -40}, {-14, -84}}, lineColor = {0, 0, 0}, textString = "mm/min")}), Documentation(info = "<html>
<p>
This block converts the input signal from metre per second to kilometre per hour and returns
the result as output signal.
</p>
</html>"));
    end from_radPerSec_to_mmPerMin;

    block to_mm
      extends Modelica.Blocks.Interfaces.PartialConversionBlock(u(unit = "m"), y(unit = "mm"));
    equation
      y = 1000 * u;
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Text(extent = {{0, 82}, {-96, 42}}, lineColor = {0, 0, 0}, textString = "m"), Text(extent = {{92, -40}, {-14, -84}}, lineColor = {0, 0, 0}, textString = "mm")}), Documentation(info = "<html>
<p>
This block converts the input signal from metre per second to kilometre per hour and returns
the result as output signal.
</p>
</html>"));
    end to_mm;
  end HelpBlocks;

  annotation(uses(
      Modelica(version="3.2.1"),
      Modelica_LinearSystems2(version="2.3.2")));
end FeedDriveLibrary;
