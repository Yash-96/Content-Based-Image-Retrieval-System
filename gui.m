function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 29-Jun-2017 13:08:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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

% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using gui.
if strcmp(get(hObject,'Visible'),'off')

end

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
queryImageFeature=[];
db=[];
c=1;
for p=1:10
    for b=0:9;
filename=randi([b*100 (b*100)+99]);
    randStrName = int2str(filename);
    randStrName = strcat(randStrName, '.jpg');
    [pathstr, name, ext] = fileparts(randStrName); % fiparts returns char type
     n=str2num(name)
     for y=0:9
         if (n>= y*100&& n<=((y*100)+99))
             db(c,y+1)=1;
         end
     end        
        a=imread(randStrName);
        handles.queryx=a;
            queryImage = imresize(a, [256 256]);
            hsvHist = hsvHistogram(queryImage);
            color_moments = colorMoments(queryImage);
            queryImage=rgb2gray(queryImage);
            lbp = extractLBPFeatures(queryImage,'Upright',true);
            glcms = graycomatrix(queryImage);
            stats = graycoprops(glcms);
            queryImageFeature(c,:) =  [color_moments hsvHist lbp stats.Contrast stats.Correlation stats.Energy stats.Homogeneity];

      c=c+1;  
    end      
end
handles.db=db;
handles.qi=queryImageFeature;
        guidata(hObject, handles);  

% --- Executes on button press in Load_Database.
function Load_Database_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=exist('database.txt');

if a>0
    msgbox('Database Already Loaded');
else
fid = fopen('database.txt', 'w+');
for i=0:999
    a=num2str(i);
    b='.jpg';
    filename=strcat(a,b);
    fprintf(fid,'%s\r',filename);
end
fclose(fid);
helpdlg('Database succesfully loaded...');
end

% --- Executes oguguin button press in Clear.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete database.txt;
% delete *.txt;
warndlg('database deleted successfully');

% --- Executes on button press in Search.
function Search_Callback(hObject, eventdata, handles)
% hObject    handle to Search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
op=handles.op;
feat=handles.feature;
db=handles.db;
queryImageFeature =handles.qi;
switch op
    case 1
    queryImageFeature=queryImageFeature(:,1:9);
    case 2
    queryImageFeature=queryImageFeature(:,10:41);
    case 3
    queryImageFeature=queryImageFeature(:,42:100);
    case 4
    queryImageFeature=queryImageFeature(:,101:104);
    case 5
    queryImageFeature=queryImageFeature(:,10:100);    
    otherwise
end  
whos('queryImageFeature');
distance=handles.dist;
queryImageFeature1=[];
d=[];    


%--------------- Open database txt file... for reading...
resultValues = zeros(1,1000);      % Results matrix..
%db1=zeros(1,1000);
%fid = fopen('database.txt');
%for z = 1:1000;
    %imagename = fgetl(fid);
    %if ~ischar(imagename), break, end       % Meaning: End of File...
     %[pathstr, name, ext] = fileparts(imagename); 
     %n1=str2double(name);
	 %db1(z)=n1;
     %X = imread(imagename);
     %queryImage1 = imresize(X, [384 256]);
        
          %  hsvHist = hsvHistogram(queryImage1);
         %   %autoCorrelogram = colorAutoCorrelogram(queryImage);
        %    color_moments = colorMoments(queryImage1);
       %     %wavelet_moments = waveletTransform(queryImage1);
      %      %img = double(rgb2gray(queryImage1))/255;
     %       %[meanAmplitude, msEnergy] = gaborWavelet(img, 4, 6);
    %        queryImage1=rgb2gray(queryImage1);
   %         lbp = extractLBPFeatures(queryImage1,'Upright',false);
  %          queryImageFeature1(z,:) = color_moments ;  
 %   end
%dlmwrite('file1.txt',db1);  
%dlmwrite('file.txt',queryImageFeature1);
%fclose(fid);
switch op
    case 1
    a=exist('color_moments1.txt');    
    if a>0
    queryImageFeature1 = dlmread('color_moments1.txt');
    end
        if a<=0
    fid = fopen('database.txt');
for z = 1:1000;
    imagename = fgetl(fid);
    if ~ischar(imagename), break, end       % Meaning: End of File...
     [pathstr, name, ext] = fileparts(imagename); 
     n1=str2double(name);
     X = imread(imagename);
     queryImage1 = imresize(X, [256 256]);
         color_moments = colorMoments(queryImage1);
         queryImageFeature1(z,:) = color_moments ;  
end
dlmwrite('color_moments1.txt',queryImageFeature1);
fclose(fid);
    end
    case 2
    a=exist('HSV1.txt');    
    if a>0
    queryImageFeature1 = dlmread('HSV1.txt');
    end
        if a<=0
    fid = fopen('database.txt');
for z = 1:1000;
    imagename = fgetl(fid);
    if ~ischar(imagename), break, end       % Meaning: End of File...
     [pathstr, name, ext] = fileparts(imagename); 
     n1=str2double(name);
     X = imread(imagename);
     queryImage1 = imresize(X, [256 256]);
        
     hsvHist = hsvHistogram(queryImage1);
          queryImageFeature1(z,:) = hsvHist ;  
   end
dlmwrite('HSV1.txt',queryImageFeature1);
fclose(fid);
    end
    case 3
    a=exist('LBP2.txt');    
    if a>0
    queryImageFeature1 = dlmread('LBP2.txt');
    end
    if a<=0
    fid = fopen('database.txt');
for z = 1:1000;
    imagename = fgetl(fid);
    if ~ischar(imagename), break, end       % Meaning: End of File...
     [pathstr, name, ext] = fileparts(imagename); 
     n1=str2double(name);
     X = imread(imagename);
     queryImage1 = imresize(X, [256 256]);
        queryImage1=rgb2gray(queryImage1);
         lbp = extractLBPFeatures(queryImage1,'Upright',true);
          queryImageFeature1(z,:) = lbp ;  
   end
dlmwrite('LBP2.txt',queryImageFeature1);
fclose(fid);
    end
    case 4
    a=exist('glcm.txt');    
    if a>0
    queryImageFeature1 = dlmread('glcm.txt');
    end
        if a<=0
    fid = fopen('database.txt');
for z = 1:1000;
    imagename = fgetl(fid);
    if ~ischar(imagename), break, end       % Meaning: End of File...
     [pathstr, name, ext] = fileparts(imagename); 
     n1=str2double(name);
     X = imread(imagename);
     queryImage1 = imresize(X, [384 256]);
                 queryImage1=rgb2gray(queryImage1);
                 glcms = graycomatrix(queryImage1);
            stats = graycoprops(glcms);
         queryImageFeature1(z,:) = [stats.Contrast stats.Correlation stats.Energy stats.Homogeneity] ;  
end
dlmwrite('glcm.txt',queryImageFeature1);
fclose(fid);
    end
    case 5
    a=exist('h.txt');    
    if a>0
    queryImageFeature1 = dlmread('h.txt');
    end
    if a<=0
    fid = fopen('database.txt');
for z = 1:1000;
    imagename = fgetl(fid);
    if ~ischar(imagename), break, end       % Meaning: End of File...
     [pathstr, name, ext] = fileparts(imagename); 
     n1=str2double(name);
     X = imread(imagename);
     queryImage1 = imresize(X, [256 256]);
           hsvHist = hsvHistogram(queryImage1);
        queryImage1=rgb2gray(queryImage1);
         lbp = extractLBPFeatures(queryImage1,'Upright',true);
          queryImageFeature1(z,:) = [hsvHist lbp];
          
end
dlmwrite('h.txt',queryImageFeature1);
fclose(fid);
    end
end
whos('queryImageFeature1');
db1=dlmread('file1.txt');
pre=zeros(1,50);
re=zeros(1,50);    
f=zeros(10,10);
for p=1:100
			for z=1:1000;
            switch distance
                case 1
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'euclidean');
                case 2
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'cityblock');
                case 3
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'chebychev');
                case 4
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'minkowski',3);
                case 5
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'cosine');
                case 6
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'correlation');
                case 7
                resultValues(z) = pdist2(queryImageFeature(p,:),queryImageFeature1(z,:),'spearman');
            end

			end
[sortedValues1, index1] = sort(resultValues);
for x=1:50
db2=zeros(1,10);
for y = 1:x;                                   
    n=db1(index1(y));
    
            for o=0:9
                if (n>= o*100&& n<=((o*100)+99))
                db2(1,o+1)= db2(1,o+1) + 1;
                end
            end
        
        
end
d(p,:)=db(p,:).*db2(1,:);
if x==20
m=floor((p-1)/10)+1;
f(m,:)=f(m,:)+d(p,:);
end
    
pre(x)=pre(x)+(sum(d(p,:))/(x));
re(x)=re(x)+(sum(d(p,:))/100);
end


end
disp(pre(20));
f1=sum(f);
f1=f1/2;
s=feat(op);
str=strjoin(s);
if op==5
figure(3);
hold on;
bar(f1);
z = mean (f1);
xl=get(gca,'Xlim');
L1=line(xl,[z z]);%mean line
set(L1,'Color','g','LineWidth',2)
title('Precision for each Class');
xlabel('Classes'); % x-axis label
ylabel('Precision'); % y-axis label
end
figure(1);
hold all;
plot(pre/100,'LineWidth',2,'DisplayName',str);
legend('-DynamicLegend');
title('Precision');
xlabel('No. of retrieved image'); % x-axis label
ylabel('Precision');
figure(2);
hold all;
plot(re/100,'LineWidth',2,'DisplayName',str);
legend('-DynamicLegend');
title('Recall');
xlabel('No. of retrieved image'); % x-axis label
ylabel('Recall');



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
distance=get(hObject,'Value');
handles.dist=distance;
contents = cellstr(get(hObject,'String'));
handles.contents=contents;
guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
distance=get(hObject,'Value');
handles.op=distance;
contents = cellstr(get(hObject,'String'));
handles.feature=contents;
guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
