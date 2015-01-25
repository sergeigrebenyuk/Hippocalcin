
function varargout = fitgui(varargin)
    %FITGUI M-file for fitgui.fig
    %      FITGUI, by itself, creates a new FITGUI or raises the existing
    %      singleton*.
    %
    %      H = FITGUI returns the handle to a new FITGUI or the handle to
    %      the existing singleton*.
    %
    %      FITGUI('Property','Value',...) creates a new FITGUI using the
    %      given property value pairs. Unrecognized properties are passed via
    %      varargin to fitgui_OpeningFcn.  This calling syntax produces a
    %      warning when there is an existing singleton*.
    %
    %      FITGUI('CALLBACK') and FITGUI('CALLBACK',hObject,...) call the
    %      local function named CALLBACK in FITGUI.M with the given input
    %      arguments.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help fitgui

    % Last Modified by GUIDE v2.5 06-Jun-2013 18:07:51

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @fitgui_OpeningFcn, ...
                       'gui_OutputFcn',  @fitgui_OutputFcn, ...
                       'gui_LayoutFcn',  [], ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
       gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT

end
%% --- Executes just before fitgui is made visible.
function fitgui_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   unrecognized PropertyName/PropertyValue pairs from the
    %            command line (see VARARGIN)

    % Choose default command line output for fitgui
    handles.output = hObject;

    global simD1 fittedD1   % Diffusion coefficient (for simulation and fitted) for cytosol
    global simD2 fittedD2   % Diffusion coefficient (for simulation and fitted) for plasma membrane
    global simK1 fittedK1   % membrane binding rate constant
    global simK2 fittedK2   % membrane unbinding rate constant
    global noise1  noise2
    global roiSize;         % Size of FRAP'ed ROI (microns)
    global modelLength      % Length of 1D dendrite (microns)
    global deltaX           % Size of data grid cell (microns)
    global simulationTime   % Simulation time (sec)
    global deltaT           % Time step in the time grid

    global hndl;
    hndl = handles;

% Toss data to controls    
    set(handles.eDistance,'String', num2str(modelLength)); 
    set(handles.eXStep,'String', num2str(deltaX)); 
    set(handles.eSimTime,'String', num2str(simulationTime)); 
    set(handles.eTStep,'String', num2str(deltaT)); 
    set(handles.eROISize,'String', num2str(roiSize)); 
    
    set(handles.eD1,'String', num2str(simD1)); 
    set(handles.eD2,'String', num2str(simD2)); 
    set(handles.eK1,'String', num2str(simK1)); 
    set(handles.eK2,'String', num2str(simK2)); 
    set(handles.fD1,'String', num2str(fittedD1)); 
    set(handles.fD2,'String', num2str(fittedD2)); 
    set(handles.fK1,'String', num2str(fittedK1)); 
    set(handles.fK2,'String', num2str(fittedK2)); 
    set(handles.eNoise1,'String', num2str(noise1)); 
    set(handles.eNoise2,'String', num2str(noise2)); 
% Clear all plots    
    cla(handles.plot_surf_cytosol)
    cla(handles.plot_2D_cytosol)
    cla(handles.plot_surf_membrane)
    cla(handles.plot_2D_membrane)
    
% Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fitgui wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

%% --- Outputs from this function are returned to the command line.
function varargout = fitgui_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

%% --- Executes on button press in openXLS.
function openXLS_Callback(hObject, eventdata, handles)
    % hObject    handle to openXLS (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    global X T rawT rawU1 rawU2 modelLength deltaX 

% Reading filenames for cytoplasm and membrane FRAP data

    [FileNameCyt,PathName,FilterIndex] = uigetfile('*.*','Import Cytoplasm FRAP Data','cyt.csv');
    [FileNameMem,PathName,FilterIndex] = uigetfile('*.*','Import Membrane FRAP Data','mem.csv');
    
    cyt = importdata(cat(2,PathName,FileNameCyt));
    mem = importdata(cat(2,PathName,FileNameMem));

% We assume same time grid for cyt and mem data        

    rawT=[];       rawT = cyt(:,1);
    set(handles.eTSize,'String', num2str(length(T))); 
    T=rawT;
% Create spatial grid

    X=[];
    for i= 1 : modelLength/ deltaX 
        X(i)= i*deltaX - modelLength/2;
    end

% Assign data    

    rawU1=[];   rawU1 = cyt(:,2);
    rawU2=[];   rawU2 = mem(:,2);

% Draw data

    DrawData(handles, 'Loaded');

end

function DrawData(handles,state)
    global Surf1  fittedSurf1
    global Surf2  fittedSurf2
    global simU1  fittedU1 rawU1 fakedU1
    global simU2  fittedU2 rawU2 fakedU2
    global X T simT rawT
    et = 0;
	switch state
                    
        case 'simulation'   % Plot simulation results
            surf(handles.plot_surf_cytosol, X,T, Surf1,'edgecolor','none')
            surf(handles.plot_surf_membrane, X,T, Surf2,'edgecolor','none')
           
            et = length(simU1); 
        
        case 'fitting'      % Plot fitting results
            surf(handles.plot_surf_cytosol, X,T, fittedSurf1,'edgecolor','none')
            surf(handles.plot_surf_membrane, X,T, fittedSurf2,'edgecolor','none')
           
            et = length(rawU1);
        
        case 'loaded'      
           et = length(rawU1);
        
        case {'one','two'}
        case 'clear'
            cla(handles.plot_surf_cytosol)
            cla(handles.plot_2D_cytosol)
            cla(handles.plot_surf_membrane)
            cla(handles.plot_2D_membrane)        
    end
    
    cla(handles.plot_2D_cytosol)
    cla(handles.plot_2D_membrane)        
    
    %hold(handles.plot_2D_cytosol)
    %hold(handles.plot_2D_membrane)
        
        if length(rawU1)% == et
            plot(handles.plot_2D_cytosol,  rawT, rawU1, '-', 'Color',[0,0,0.7], 'LineWidth',2)
            plot(handles.plot_2D_membrane, rawT, rawU2, '-', 'Color',[0,0,0.7], 'LineWidth',2)
        end
        if length(simU1) %== et
            plot(handles.plot_2D_cytosol,  simT, simU1, '-', 'Color',[0,0.5,0.5], 'LineWidth',2)
            plot(handles.plot_2D_membrane, simT, simU2, '-', 'Color',[0,0.5,0.5], 'LineWidth',2)    
        end
        if length(fakedU1)% == et
            plot(handles.plot_2D_cytosol,  simT, fakedU1,'.', 'Color',[0,0.5,0.5],'LineWidth',2, 'MarkerSize',10)
            plot(handles.plot_2D_membrane, simT, fakedU2,'.', 'Color',[0,0.5,0.5],'LineWidth',2, 'MarkerSize',10)
        end
        if length(fittedU1)% == et
            plot(handles.plot_2D_cytosol,  rawT, fittedU1,'-', 'Color',[1,0,0], 'LineWidth',2)
            plot(handles.plot_2D_membrane, rawT, fittedU2,'-', 'Color',[1,0,0], 'LineWidth',2)
        end
   % hold off
    
end

% Not used
function stop = OptStep(x, optimValues, state)
    global IterCnt;
    IterCnt = IterCnt +1;
    global hndl;
    
    set(hndl.eIteration,'String', num2str(IterCnt)); 
    guidata(hndl.eIteration, hndl);
    stop = false;
end
function [u] = fitfun(v,xdata)
    global simD2;
    [s1,s2,u1,u2]=pdesolve(v(1),v(2),simD2,v(3));
    u = [u1;u2];
end
%% --- Executes on button press in DoFit.
function DoFit_Callback(hObject, eventdata, handles)
% hObject    handle to DoFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global fittedD1   % Diffusion coefficient (for simulation and fitted) for cytosol
    global fittedD2   % Diffusion coefficient (for simulation and fitted) for plasma membrane
    global fittedK1   % membrane binding rate constant
    global fittedK2   % membrane unbinding rate constant
    global Surf1  fittedSurf1
    global Surf2  fittedSurf2
    global simU1  fittedU1 rawU1 fakedU1
    global simU2  fittedU2 rawU2 fakedU2
    global IterCnt;
    global roiSize;         % Size of FRAP'ed ROI (microns)
    global modelLength      % Length of 1D dendrite (microns)
    global deltaX           % Size of data grid cell (microns)
    global simulationTime   % Simulation time (sec)
    global deltaT           % Time step in the time grid
    global X   % Vector [x0, x1, ..., xn] specifying the points at which a numerical solution is requested for every value in tspan ds
    global T   % Vector [t0, t1,..., tf] specifying the points at which a solution is requested for every value in X
    
    set(handles.plot_2D_cytosol,'Color', [1,0.9,0.9]); set(handles.plot_2D_membrane,'Color', [1,0.9,0.9]); 

%% Read stating values of parameters
    
    fittedD1 = str2double(get(handles.fD1, 'string'));
    fittedD2 = str2double(get(handles.fD2, 'string'));
    fittedK1 = str2double(get(handles.fK1, 'string'));
    fittedK2 = str2double(get(handles.fK2, 'string'));
    modelLength = str2double(get(handles.eDistance, 'string'));
    deltaX = str2double(get(handles.eXStep, 'string'));

%%  Checking solution grid
    
	if (~length(T)) 
        return;
    end; 
    X=[];
    for i= 1 : modelLength/ deltaX 
        X(i)= i*deltaX - modelLength/2;
        %X(i)= i*deltaX;
    end


%%  Fitting our model to data

%   Start point [D1, k1, k2]
    v0 = [fittedD1,fittedK1,fittedK2];
%        [D1  k1  k2 ]
    lb = [0,  0,  0  ];
    ub = [200,100,100];
    IterCnt = 0;
    options = optimoptions(@lsqcurvefit,'PlotFcns',{@optimplotx;@optimplotresnorm;@optimplotstepsize},...
        'TolX',1e-06 ,...
        'FinDiffType','central');
      %  'OutputFcn',@OptStep);
    %options = optimoptions(@lsqcurvefit,'OutputFcn',@OptStep);
    raw_data = [rawU1;rawU2];
    [x,resnorm] = lsqcurvefit(@fitfun,v0,T,raw_data,lb,ub,options);
    
    fittedD1 = x(1);
    fittedK1 = x(2);
    fittedK2 = x(3);
    

%% Generate report
    
    set(handles.fD1,'String', num2str(fittedD1)); 
    set(handles.fK1,'String', num2str(fittedK1)); 
    set(handles.fK2,'String', num2str(fittedK2)); 
    
    [fittedSurf1,fittedSurf2,fittedU1,fittedU2]=pdesolve(x(1),x(2),fittedD2,x(3));
    set(handles.plot_2D_cytosol,'Color', [1,1,1]); set(handles.plot_2D_membrane,'Color', [1,1,1]); 
    DrawData(handles,'fitting');
end

% --- Executes on button press in bSimulate.
function bSimulate_Callback(hObject, eventdata, handles)
% hObject    handle to bSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global simD1 % Diffusion coefficient (for simulation and fitted) for cytosol
    global simD2 % Diffusion coefficient (for simulation and fitted) for plasma membrane
    global simK1 % membrane binding rate constant
    global simK2 % membrane unbinding rate constant
    global Surf1
    global Surf2
    global simU1 fakedU1
    global simU2 fakedU2
    global noise1  noise2
    global roiSize;         % Size of FRAP'ed ROI (microns)
    global modelLength      % Length of 1D dendrite (microns)
    global deltaX           % Size of data grid cell (microns)
    global simulationTime   % Simulation time (sec)
    global deltaT           % Time step in the time grid
    global X T simT   % Vectors [t0, t1,..., tf] specifying the points at which a solution is requested for every value in X

% Read user parameters from edit boxes

    simD1 = str2double(get(handles.eD1, 'string'));
    simD2 = str2double(get(handles.eD2, 'string'));
    simK1 = str2double(get(handles.eK1, 'string'));
    simK2 = str2double(get(handles.eK2, 'string'));
    noise1 = str2double(get(handles.eNoise1, 'string'));
    noise2 = str2double(get(handles.eNoise2, 'string'));
    
    roiSize = str2double(get(handles.eROISize, 'string'));
    modelLength = str2double(get(handles.eDistance, 'string'));
    deltaX = str2double(get(handles.eXStep, 'string'));
    simulationTime = str2double(get(handles.eSimTime, 'string'));
    deltaT = str2double(get(handles.eTStep, 'string'));

% Prepare solution grid

    X=[];
    for i= 1 : modelLength/ deltaX 
        X(i)= i*deltaX - modelLength/2;
        %X(i)= i*deltaX;
    end
    simT=[];
    for i=1 : simulationTime / deltaT
        simT(i)=i*deltaT;
    end
    T=simT;
% Simulate data : simU1, simU2 will contain solution at x=0;

    [Surf1,Surf2,simU1,simU2]=pdesolve(simD1,simK1,simD2,simK2);

% Add some noise to simulated data

    fakedU1 = simU1 + noise1*randn(simulationTime / deltaT,1);
    fakedU2 = simU2 + noise2*randn(simulationTime / deltaT,1);

% Draw data
 
    DrawData(handles,'simulation');
%    surf(handles.plot_surf_cytosol, X,T, Surf1,'edgecolor','none')
%    plot(handles.plot_2D_cytosol, T, simU1, '-', T, fakedU1,'.','LineWidth',2,'MarkerSize',10);
%    surf(handles.plot_surf_membrane, X,T, Surf2,'edgecolor','none')
%    plot(handles.plot_2D_membrane, T, simU2, '-',T, fakedU2,'.','LineWidth',2,'MarkerSize',10);
   
end


% --- Executes on button press in bMakeUpData.
function bMakeUpData_Callback(hObject, eventdata, handles)
% hObject    handle to bMakeUpData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global rawU1 fakedU1
    global rawU2 fakedU2
    global T rawT
    rawU1 = fakedU1;
    rawU2 = fakedU2;

    T = rawT;
end


function eD1_Callback(hObject, eventdata, handles)
% hObject    handle to eD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eD1 as text
%        str2double(get(hObject,'String')) returns contents of eD1 as a double
end

% --- Executes during object creation, after setting all properties.
function eD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function eK1_Callback(hObject, eventdata, handles)
% hObject    handle to eK1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eK1 as text
%        str2double(get(hObject,'String')) returns contents of eK1 as a double

end
% --- Executes during object creation, after setting all properties.
function eK1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eK1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function eD2_Callback(hObject, eventdata, handles)
% hObject    handle to eD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eD2 as text
%        str2double(get(hObject,'String')) returns contents of eD2 as a double
end

% --- Executes during object creation, after setting all properties.
function eD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function eK2_Callback(hObject, eventdata, handles)
% hObject    handle to eK2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eK2 as text
%        str2double(get(hObject,'String')) returns contents of eK2 as a double

end
% --- Executes during object creation, after setting all properties.
function eK2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eK2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function eDistance_Callback(hObject, eventdata, handles)
% hObject    handle to eDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eDistance as text
%        str2double(get(hObject,'String')) returns contents of eDistance as a double
end

% --- Executes during object creation, after setting all properties.
function eDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function eXStep_Callback(hObject, eventdata, handles)
% hObject    handle to eXStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eXStep as text
%        str2double(get(hObject,'String')) returns contents of eXStep as a double
end

% --- Executes during object creation, after setting all properties.
function eXStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eXStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function eSimTime_Callback(hObject, eventdata, handles)
% hObject    handle to eSimTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eSimTime as text
%        str2double(get(hObject,'String')) returns contents of eSimTime as a double
end

% --- Executes during object creation, after setting all properties.
function eSimTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eSimTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function eTStep_Callback(hObject, eventdata, handles)
% hObject    handle to eTStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eTStep as text
%        str2double(get(hObject,'String')) returns contents of eTStep as a double
end

% --- Executes during object creation, after setting all properties.
function eTStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eTStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function eROISize_Callback(hObject, eventdata, handles)
% hObject    handle to eROISize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eROISize as text
%        str2double(get(hObject,'String')) returns contents of eROISize as a double
end

% --- Executes during object creation, after setting all properties.
function eROISize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eROISize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function eNoise1_Callback(hObject, eventdata, handles)
% hObject    handle to eNoise1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eNoise1 as text
%        str2double(get(hObject,'String')) returns contents of eNoise1 as a double
end

% --- Executes during object creation, after setting all properties.
function eNoise1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eNoise1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function eNoise2_Callback(hObject, eventdata, handles)
% hObject    handle to eNoise2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eNoise2 as text
%        str2double(get(hObject,'String')) returns contents of eNoise2 as a double
end

% --- Executes during object creation, after setting all properties.
function eNoise2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eNoise2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes when main_fig is resized.
function main_fig_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to main_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function fD2_Callback(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text16 as text
%        str2double(get(hObject,'String')) returns contents of text16 as a double
end


% --- Executes during object creation, after setting all properties.
function fD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in bClearData.
function bClearData_Callback(hObject, eventdata, handles)
% hObject    handle to bClearData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    DrawData(handles,'clear');
end
