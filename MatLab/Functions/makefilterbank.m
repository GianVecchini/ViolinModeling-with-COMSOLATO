function [sim_values, sim_time] = makefilterbank(centralFreq, bandWidth, peakMagnitude, N, num_samples, fs, input_s, input_t)

%% filterbank constructor
close_system("Filterbank", 0)


orygin = [100 100];
xOffset = 200;
yOffset = 100;
blockSize = [50 50];

ResLbl = cell(1,N);
IndLbl = cell(1,N);
CapLbl = cell(1,N);

open_system("Filterbank");

%insert N block-lines
for n = 1:N

    %center = centralFreq(n) * peakMagnitude(n);
    center = centralFreq(n);

    syms x y z

    % x = R, y = L, z = C
    eq1 = 1/(2*pi*sqrt(y*z)) == centralFreq(n); % Central frequency
    eq2 = x/(2*pi*y) == bandWidth(n);           % BandWidth
    eq3 = 1/x == peakMagnitude(n);              % Magnitude

    solutions = solve([eq1, eq2, eq3], [x, y, z]);

    R = abs(double(solutions.x(1)));
    L = abs(double(solutions.y(1)));
    C = abs(double(solutions.z(1)));

    % set labels
    ResLbl{n} = ['R' num2str(n)];
    IndLbl{n} = ['L' num2str(n)]; 
    CapLbl{n} = ['C' num2str(n)];
    
    % set position
    resPos = [orygin(1) orygin(2) + n*yOffset orygin(1) + blockSize(1) orygin(2) + n*yOffset + blockSize(2)];
    indPos = [orygin(1) + yOffset orygin(2) + n*yOffset orygin(1) + yOffset + blockSize(1) orygin(2) + n*yOffset + blockSize(2)];
    capPos = [orygin(1) + 2*yOffset orygin(2) + n*yOffset orygin(1) + 2*yOffset + blockSize(1) orygin(2) + n*yOffset + blockSize(2)];

    % insert resistor
    blockName = ['Filterbank/' ResLbl{n}];
    add_block('fl_lib/Electrical/Electrical Elements/Resistor', blockName);
    set_param(blockName, 'Position', resPos);
    set_param(blockName,'R',num2str(R(1)));
    set_param(blockName,'AttributesFormatString', num2str(R(1)));

    % insert inductor
    blockName = ['Filterbank/' IndLbl{n}];
    add_block('fl_lib/Electrical/Electrical Elements/Inductor', blockName);
    set_param(blockName, 'Position', indPos);
    set_param(blockName,'l',num2str(L(1)));
    set_param(blockName,'AttributesFormatString', num2str(L(1)));
    

    %insert capacitor
    blockName = ['Filterbank/' CapLbl{n}];
    add_block('fl_lib/Electrical/Electrical Elements/Capacitor', blockName);
    set_param(blockName, 'Position', capPos);
    set_param(blockName,'c',num2str(C(1)))
    set_param(blockName,'r',num2str(0))
    set_param(blockName,'AttributesFormatString', num2str(C(1)));

    % serial connections
    add_line('Filterbank', [ResLbl{n} '/RConn1'], [IndLbl{n} '/LConn1']);
    add_line('Filterbank', [IndLbl{n} '/RConn1'], [CapLbl{n} '/LConn1']);
    disp(['f0 = ', num2str(centralFreq(n))]);
    disp(['R = ', num2str(R(1))]);
    disp(['L = ', num2str(L(1))]);
    disp(['C = ', num2str(C(1))]);
    disp('-------o--------');


end

% parallel connentions
for n = 1:N-1
    add_line('Filterbank', [ResLbl{n} '/LConn1'], [ResLbl{n + 1} '/LConn1'] );
    add_line('Filterbank', [CapLbl{n} '/RConn1'], [CapLbl{n + 1} '/RConn1'] )
end

%sensor
%sensPos = [capPos(1) + 200 capPos(2) capPos(3) + 200 capPos(4)];
sensPos = capPos + [200 50 200 50];
add_block('fl_lib/Electrical/Electrical Sensors/Current Sensor','Filterbank/Sensor');
set_param('Filterbank/Sensor', 'Position', sensPos);
set_param('Filterbank/Sensor', 'Orientation', 'right');
set_param('Filterbank/Sensor', 'BlockMirror', 'off'); 
add_line('Filterbank', [CapLbl{end} '/RConn1'], 'Sensor/LConn1');

% source
sourcePos = resPos + [-200 50 -200 50];
add_block('fl_lib/Electrical/Electrical Sources/Controlled Voltage Source','Filterbank/Source');
set_param('Filterbank/Source', 'Position', sourcePos);
add_line('Filterbank', 'Source/LConn1', [ResLbl{end} '/LConn1']);

%ground
groundPos = sourcePos + [300 200 300 200];
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference','Filterbank/Ground');
set_param('Filterbank/Ground', 'Position', groundPos);
add_line('Filterbank', 'Ground/LConn1', 'Source/RConn2');
add_line('Filterbank', 'Ground/LConn1', 'Sensor/RConn2');

%output
add_block('nesl_utility/PS-Simulink Converter', 'Filterbank/OutConv');
set_param('Filterbank/OutConv', 'Position', sensPos + [150 0 150 0]);
add_line('Filterbank', 'Sensor/Rconn1', 'OutConv/Lconn1');
add_block('simulink/Sinks/To Workspace', 'Filterbank/Simout');
set_param('Filterbank/Simout', 'SaveFormat', 'StructureWithTime');
set_param('Filterbank/Simout', 'SampleTime', '-1');
set_param('Filterbank/Simout', 'Position', sensPos + [250 0 250 0]);
add_line('Filterbank', 'OutConv/1', 'Simout/1');

%input
add_block('nesl_utility/Simulink-PS Converter', 'Filterbank/InConv');
set_param('Filterbank/InConv', 'Position', sourcePos + [-100 50 -100 50]);
add_line('Filterbank', 'Source/Rconn1', 'InConv/Rconn1');
add_block('simulink/Sources/From Workspace', 'Filterbank/Simin');
set_param('Filterbank/Simin', 'SampleTime', num2str(1/fs));
set_param('Filterbank/Simin', 'Position', sourcePos + [-200 50 -200 50]);
add_line('Filterbank', 'Simin/1', 'InConv/1');

%solver config
add_block('nesl_utility/Solver Configuration', 'Filterbank/Config');
set_param('Filterbank/Config', 'Position', sourcePos + [-100 -100 -100 -100]);
add_line('Filterbank', 'Source/Lconn1', 'Config/Rconn1');

%% simulate
simin = timeseries(input_s, input_t, 'Name', 'ImpulsiveIn');
assignin("base", 'simin', simin);

%% out
set_param('Filterbank', 'StopTime', num2str(input_t(end)));
simOut = sim('Filterbank');
simout = simOut.get('simout');
sim_values = simout.signals.values;
sim_time = simout.time;
end