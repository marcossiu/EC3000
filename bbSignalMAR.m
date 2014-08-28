% This function was developed for attendance of the course
% EC4530 - Software Radio
% Assignment: LAB1, section 4.11.
% written by: Marcos Siu (msiu@nps.edu)
% version 2: Aug04,2014
%
% bbSignalMAR function:
%
% This functions returns a vector that represents the samples of a complex
% envelope of baseband Signal.
%
% function bbsignal = bbSignalMAR(a,B,pulseParam,M)
% a = signal constellation values;
% E = average symbol energy;
% M = oversampling factor (2,4 or 8);
% pulseParam = Structure of variables that gives the pulse characteristics:
% Example:  pulseParam.symInterval sets the pulse length (in seconds)
%           pulseParam.type = 'rectangular' sets the pulse as a rectangular
%           pulse. pulseParam.type is a string variable.
%
% Example: Given a NRZ signal where:
%                
% a = {+1 or -1}
%                   +N
%                  _____
%                  \       
% bbSignalMAR = B*  \    a[k].p(t-kT),     0 <= k <= N-1 (defined by user)
%                   /                   
%                  /____ 
%                   k=0
%
% For a unit-energy rectangular pulse B = 1/T.
%
% v2: improvements: 
% - B (voltage level in function of E (pulse energy) (assignment 2.1);
% - D and alpha reading (duration of the pulse and rolloff);
% - Buffers I and Q are expanded:
% in order to support first from 0 to 0.5DT and 
% last symbol from (N-1)T + 0.5DT.


function bbsignal = bbSignalMAR(a,E,pulseParam,M)
% EC4530: Software Radio
% Students: Marcos Siu - msiu@nps.edu
%
% LAB1) SDR transmitter that modulates the signal space points.
% Section: 4.12 Code
% PART1) bbSignal Function
%
%Variables:
T = pulseParam.symInterval; %pulse interval (period)
type = pulseParam.type; %Pulse type (rectangular, etc.)
alpha = pulseParam.rolloff; %Roll Off factor
D = pulseParam.durInSym; %the duration of the pulse in symbols (must be even value).

% Check parameters integrity:
if T < 0 | rem(M,1) > eps | rem(D,2) > eps
         %display an error return:
        disp('arguments mismatch')
        return
end

Tsamp = T/M; %Sample period (<Sample Frequency>^(-1))
N = length(a)+D; %computes the MSG size.
n = 0:1:((M*N)-1); %sample size (M*N)
t = n*Tsamp; %time defined in function of Sample period or Sample frequency.
buffer_I = zeros(N-D,N*M); %Temp buffer to store pulses for in-phase pulses.
buffer_Q = zeros(N-D,N*M); %Temp buffer to store pulses for quadrature pulses.
k = 0; %sum index
B = sqrt(E); %equivalent voltage to achieve the desired pulse energy.

for k = 0:1:(N-D-1)
    %the following step stores the k-th pulse in the k-th row.
    buffer_I(k+1,:) = pulseMAR(t-(k*T),pulseParam);
    %Stores symbol amplitude in-phase buffer
    buffer_I(k+1,:) = real(a(1,k+1)).* buffer_I(k+1,:); 
    if imag(a(1,k+1)) ~= 0;
        %the following step stores the k-th pulse in the k-th row.
        buffer_Q(k+1,:) = pulseMAR(t-(k*T),pulseParam);
        %Stores symbol amplitude in quadrature buffer
        buffer_Q(k+1,:) = imag(a(1,k+1)).* buffer_Q(k+1,:);
    end
end

%combining in-phase and quadrature signal
bbsignal = B*sum(buffer_I) + (1j*B)*sum(buffer_Q);

return
    