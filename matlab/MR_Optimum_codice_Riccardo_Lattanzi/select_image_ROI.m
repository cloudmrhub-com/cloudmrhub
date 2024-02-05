function varargout = select_image_ROI(varargin)
% SELECT_IMAGE_ROI M-file for select_image_ROI.fig
%      SELECT_IMAGE_ROI, by itself, creates a new SELECT_IMAGE_ROI or raises the existing
%      singleton*.
%
%      H = SELECT_IMAGE_ROI returns the handle to a new SELECT_IMAGE_ROI or the handle to
%      the existing singleton*.
%
%      SELECT_IMAGE_ROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_IMAGE_ROI.M with the given input arguments.
%
%      SELECT_IMAGE_ROI('Property','Value',...) creates a new SELECT_IMAGE_ROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_image_ROI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_image_ROI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help select_image_ROI

% Last Modified by GUIDE v2.5 14-Mar-2005 17:14:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_image_ROI_OpeningFcn, ...
                   'gui_OutputFcn',  @select_image_ROI_OutputFcn, ...
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

%--------------------------------------------------------------------
% --- Executes just before select_image_ROI is made visible.
function select_image_ROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_image_ROI (see VARARGIN)

% Choose default command line output for select_image_ROI
handles.output = [];
guidata(hObject, handles);
set(hObject, 'WindowStyle', 'modal');

axes(handles.axes_image);

handles.image_data.y_dim = varargin{1};
handles.image_data.z_dim = varargin{2};
set(handles.fig_select, 'NextPlot', 'Add');
one_slice_image = varargin{3};
imshow(abs(one_slice_image), []);
guidata(hObject, handles);

% need another transparent overlaying axes because ButtonDownFcn callback
% does not execute when click over on another graphics object displayed in
% the axes
axes(handles.axes_overlay);
set(handles.axes_overlay, 'Color', 'none');
set(handles.axes_overlay, 'XTick', [], 'YTick', []);
box on;
% Lock the axis so it will not resize 
axis manual;

handles.fig_data.step_size = ResizeAxesOverlay(handles);
guidata(hObject, handles);

handles.fig_data.regions = [];
set(handles.listbox_region, 'String', {});
handles.fig_data.selected_region = 0;
handles.fig_data.rect_count = 1;
% Store a copy of the positions of the rectangles in a matrix for rapid
% processing when handling clicks and drags.
handles.fig_data.rect_pos = [];

if length(varargin) == 4
    input_param = varargin{4};
    num_rect = size(input_param);
    if num_rect(2) == 4
        % If the input matrix has four columns, draw the rectangles according
        % to the input

        % for a n-by-m matrix, imshow's orientation of the image is
        % (1,1) ---------------- (m,1) 
        %       |              |
        %       |              |
        %       |              |
        %       |              |
        % (1,n) ---------------- (m,n)

        % The line is drawn to the right (and top) side of the pixel. To include
        % the leftmost and buttom most pixels of the selection in the rectangle
        % substract one from the minimum such that the line is draw on the
        % right side and top side of the minus-1 pixels (left side and buttom
        % side of the selection)
        y_dim = handles.image_data.y_dim;
        z_dim = handles.image_data.z_dim;
        input_rect = zeros(num_rect);
        input_rect(:,1) = (input_param(:,1) - 1)./y_dim;
        input_rect(:,2) = 1 - input_param(:,4)./z_dim;
        input_rect(:,3) = (input_param(:,2) - input_param(:,1) + 1)./y_dim;
        input_rect(:,4) = (input_param(:,4) - input_param(:,3) + 1)./z_dim;



        for ctr = 1:num_rect(1)
            if ~OverlapOtherRegions(handles.fig_data.rect_pos, input_rect(ctr,:))
                handles = DrawNewRectangle(handles, input_rect(ctr,:));
            else
                set(handles.text_instruction, 'String', ...
                    'Some regions specified by input overlaps each other; those regions are not displayed.');
            end
        end
        SetSelectedRegion(handles);
    end
else
    set(handles.text_instruction, 'String', ...
    'To add a region, click on a point not within any regions and drag the mouse.');
end

guidata(hObject, handles);
% UIWAIT makes select_image_ROI wait for user response (see UIRESUME)
uiwait(handles.fig_select);


%--------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = select_image_ROI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(gcf);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Axes Panel           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
% --- Executes on mouse press over axes background.
function axes_overlay_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_overlay);
point1 = get(gca, 'CurrentPoint');
point1 = point1(1, 1:2);

% Round the points to the nearest pixel (in the matrix) so that in
% the final rectangle drawn, there will be no partial pixels
y_dim = handles.image_data.y_dim;
z_dim = handles.image_data.z_dim;
point1(1) = round(point1(1)*y_dim)/y_dim;
point1(2) = round(point1(2)*z_dim)/z_dim;

% If the click is not within the axes, do nothing
if WithinAxes(point1)
    % Since regions do not overlap each other, a point can be within
    % only one region
    index_region = WithinWhichRegion(handles.fig_data.rect_pos, point1);
    
    % If the point is not within any region, draw a new rectangle
    if isempty(index_region)
        UnsetSelectedRegion(handles);
        
        % If the click is not within a selected region, this click will 
        % produce a rubberband box for use to selecte a region
        set(handles.text_instruction, 'String', ...
            'Drag the box to select a new region');
        finalRect = rbbox;
        point2 = get(gca, 'CurrentPoint');
        point2 = point2(1, 1:2);
        point2(1) = round(point2(1)*y_dim)/y_dim;
        point2(2) = round(point2(2)*z_dim)/z_dim;
        
        start_point = min(point1, point2);
        side_dim = abs(point1 - point2);
        
        if WithinAxes(point2) && side_dim(1) > 0 && side_dim(2) > 0
            % If the rectangle is within the axes and the selected region
            % has width and length greater than 0, proceed to other
            % checks, Only need to check point2 because it is impossible 
            % for point1 to lie outside the axes (callback would not
            % be executed)
           new_rect_pos = horzcat(start_point, side_dim);
 
            % If the new region does not overlap other existing regions,
            % draw the rectangle; otherwise give a warning in the instruction
            % text and abort.
            if ~OverlapOtherRegions(handles.fig_data.rect_pos, new_rect_pos)
                handles = DrawNewRectangle(handles, new_rect_pos);
            else
                set(handles.text_instruction, 'String', ...
                    'New region must not overlap existing ones');
            end
        end
        
        SetSelectedRegion(handles);
        
    else
        % If the click is within an existing region, highlight it and
        % display a dragrect box
        if index_region ~= handles.fig_data.selected_region
            UnsetSelectedRegion(handles);
            handles.fig_data.selected_region = index_region;
            SetSelectedRegion(handles);
        end
        
        set(handles.text_instruction, 'String', ...
            'Drag the region to a new position.');
        % Convert the position vector from normalized to pixel (for
        % dragrect)
        region_pos = get(handles.fig_data.regions(index_region), 'Position');
        axes_pos = get(handles.axes_overlay, 'Position');
        fig_pos = get(handles.fig_select, 'Position');
        % Convert the axes position vector into pixel unit
        axes_pos = [ axes_pos(1)*fig_pos(3) axes_pos(2)*fig_pos(4) ...
            axes_pos(3)*fig_pos(3) axes_pos(4)*fig_pos(4) ];
        % Conver the enclosing rectangle position vector into pixel unit
        region_pos_pixel = [ region_pos(1)*axes_pos(3) + axes_pos(1) ...
            region_pos(2)*axes_pos(4) + axes_pos(2) ...
            region_pos(3)*axes_pos(3) region_pos(4)*axes_pos(4) ];
        [ new_region_pos ] = dragrect(region_pos_pixel);

        % Convert the pixels in new_region_pos back into normalized unit
        % (for drawing the rectangle in axes)
        new_region_pos = [ (new_region_pos(1) - axes_pos(1))/axes_pos(3) ...
            (new_region_pos(2) - axes_pos(2))/axes_pos(4) ...
            region_pos(3) region_pos(4) ];
        % Snap it to the pixel grid
        new_region_pos(1) = round(new_region_pos(1)*y_dim)/y_dim;
        new_region_pos(2) = round(new_region_pos(2)*z_dim)/z_dim;
        % Only move the rectangle if it is within the limits of the axes
        % and does not overlap other regions (minus itself)
        rect_pos = handles.fig_data.rect_pos;
        rect_pos(index_region, :) = [];
        if new_region_pos(1) >= 0 && new_region_pos(2) >= 0 && ...
                new_region_pos(1) + new_region_pos(3) <= 1 && ...
                new_region_pos(2) + new_region_pos(4) <= 1 && ...
                ~OverlapOtherRegions(rect_pos, new_region_pos)
            set(handles.fig_data.regions(index_region), 'Position', new_region_pos);
            handles.fig_data.rect_pos(index_region,:) = new_region_pos;
            set(handles.text_instruction, 'String', '');
        else
            set(handles.text_instruction, 'String', ...
                'Region must stay within the image and does not overlap with other regions.');
        end
    end
end

guidata(hObject, handles);




%--------------------------------------------------------------------
function bool = WithinAxes(point)
bool = (point(1) >= 0 && point(1) <= 1 && point(2) >= 0 && point(2) <= 1);

%--------------------------------------------------------------------
function index_region = WithinWhichRegion(rect_pos, point)
if isempty(rect_pos)
    index_region = [];
else
    num_rect = size(rect_pos);
    % Translate the rect_pos into [ y_min y_max z_min z_max; ... ]
    rect_pos = [ rect_pos(:,1) rect_pos(:,1) + rect_pos(:,3) ...
        rect_pos(:,2) rect_pos(:,2) + rect_pos(:,4) ];

    % Create a matrix the same size as rect_pos and compare the columns
    mat_compare = repmat([ point(1) point(1) point(2) point(2) ], [ num_rect(1) 1 ]);
    mat_compare(:,1) = ( mat_compare(:,1) > rect_pos(:,1) );
    mat_compare(:,2) = ( mat_compare(:,2) < rect_pos(:,2) );
    mat_compare(:,3) = ( mat_compare(:,3) > rect_pos(:,3) );
    mat_compare(:,4) = ( mat_compare(:,4) < rect_pos(:,4) );
    index_region = find(all(mat_compare, 2));
end

%--------------------------------------------------------------------
function bool = OverlapOtherRegions(rect_pos, new_rect)
bool = 1;
if isempty(rect_pos)
    bool = 0;
else
    % For two rectangles (x1, y1, h1, w1) and (x2, y2, h2, w2), they do not
    % overlap each other if 
    %   x1 + h1 < x2    or
    %   x2 + h2 < x1    or
    %   y1 + w1 < y2    or
    %   y2 + w2 < y1
    % if none of these holds then they overlap
    
    % Calculate the matrix where [ x1+h1 x2+h2 y1+w1 y2+w2; ... ]
    mat_compare = zeros(size(rect_pos));
    mat_compare(:,1) = rect_pos(:,1) + rect_pos(:,3);
    mat_compare(:,2) = new_rect(1) + new_rect(3);
    mat_compare(:,3) = rect_pos(:,2) + rect_pos(:,4);
    mat_compare(:,4) = new_rect(2) + new_rect(4);
    % Compare with [ x2 x1 y2 y1 ]
    mat_compare(:,1) = ( mat_compare(:,1) <= new_rect(1) );
    mat_compare(:,2) = ( mat_compare(:,2) <= rect_pos(:,1) );
    mat_compare(:,3) = ( mat_compare(:,3) <= new_rect(2) );
    mat_compare(:,4) = ( mat_compare(:,4) <= rect_pos(:,2) );
    
    if all(any(mat_compare, 2))
        bool = 0;
    end
end

%--------------------------------------------------------------------
function handles = DrawNewRectangle(handles, rect_pos)
len = length(handles.fig_data.regions);
handles.fig_data.regions(len+1) = rectangle('Position', rect_pos, ...
    'Curvature', [0, 0], 'EdgeColor', 'g');
str_listbox = get(handles.listbox_region, 'String');

% Update the listbox and the rect_pos matrix
str_listbox{len+1} = [ 'Rectangle #' num2str(handles.fig_data.rect_count) ];
handles.fig_data.rect_count = handles.fig_data.rect_count + 1;
set(handles.listbox_region, 'String', str_listbox);
set(handles.listbox_region, 'Value', len+1);
handles.fig_data.selected_region = len+1;
handles.fig_data.rect_pos(len+1,:) = rect_pos;
set(handles.text_instruction, 'String', '');

%--------------------------------------------------------------------
function UnsetSelectedRegion(handles)
lb_select = handles.fig_data.selected_region;
if length(handles.fig_data.regions) > 0 && lb_select ~= 0
    selected_region = handles.fig_data.regions(lb_select);
    set(selected_region, 'EdgeColor', 'g');
end

%--------------------------------------------------------------------
function SetSelectedRegion(handles)
lb_select = handles.fig_data.selected_region;
if lb_select ~= 0
    selected_region = handles.fig_data.regions(lb_select);
    set(selected_region, 'EdgeColor', 'r');
    set(handles.listbox_region, 'Value', lb_select);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Region Panel          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
% --- Executes on selection change in listbox_region.
function listbox_region_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_region contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_region
UnsetSelectedRegion(handles);
handles.fig_data.selected_region = get(hObject, 'Value');
guidata(hObject, handles);
SetSelectedRegion(handles);


%--------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function listbox_region_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------
% --- Executes on button press in but_delete.
function but_delete_Callback(hObject, eventdata, handles)
% hObject    handle to but_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lb_select = handles.fig_data.selected_region;
% Make the rectangle invisible before deleting it (display does not refresh
% automatically after deletion)
set(handles.fig_data.regions(lb_select), 'Visible', 'off');
handles.fig_data.regions(lb_select) = [];
handles.fig_data.rect_pos(lb_select,:) = [];
str_listbox = get(handles.listbox_region, 'String');
str_listbox(lb_select) = [];
set(handles.listbox_region, 'String', str_listbox);
if length(handles.fig_data.regions) == 0
    handles.fig_data.selected_region = 0;
else
    handles.fig_data.selected_region = length(handles.fig_data.regions);
    set(handles.listbox_region, 'Value', handles.fig_data.selected_region);
end
guidata(hObject, handles);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Size Adjustment         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
% --- Executes when fig_select is resized.
function fig_select_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to fig_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fig_data.step_size = ResizeAxesOverlay(handles);
guidata(hObject, handles);



%--------------------------------------------------------------------
function new_step_size = ResizeAxesOverlay(handles)
% Adjust axes_overlay's position such that it is still fitted around the resized
% image
y_dim = handles.image_data.y_dim;
z_dim = handles.image_data.z_dim;
axes_pos = get(handles.axes_image, 'Position');
figure_pos = get(handles.fig_select, 'Position');
image_ratio = y_dim/z_dim;
axes_pixel_ratio = (axes_pos(3)*figure_pos(3))/(axes_pos(4)*figure_pos(4));

% Determine which side of the figure constraints the size of the image
% after resizing
if axes_pixel_ratio >= image_ratio
    % Figure height is the constraint and so is set to maximum.
    % Adjust width and starting x of the overlying axes to fit the
    % underlying image
    set(handles.axes_overlay, 'Position', ...
        AdjustWidth([ y_dim z_dim ], figure_pos, axes_pos));
    new_step_size = axes_pos(3)*figure_pos(3)/y_dim;
else
    % Otherwise, figure width is the constraint. Adjust the height.
    set(handles.axes_overlay, 'Position', ...
        AdjustHeight([ y_dim z_dim ], figure_pos, axes_pos));
    new_step_size = axes_pos(4)*figure_pos(4)/z_dim;
end


%--------------------------------------------------------------------
% Both functions assume that the image is in the middle of the figure 
function new_pos = AdjustWidth(image_dim, figure_pos, axes_pos)
image_ratio = image_dim(1)/image_dim(2);
axes_pixel_ratio = (axes_pos(4)*figure_pos(4))/(axes_pos(3)*figure_pos(3));
new_width = axes_pos(3) * axes_pixel_ratio * image_ratio;
new_x = axes_pos(1) + (axes_pos(3) - new_width)/2;
new_pos = [ new_x axes_pos(2) new_width axes_pos(4) ];


%--------------------------------------------------------------------
function new_pos = AdjustHeight(image_dim, figure_pos, axes_pos)
image_ratio = image_dim(2)/image_dim(1);
axes_pixel_ratio = (axes_pos(3)*figure_pos(3))/(axes_pos(4)*figure_pos(4));
new_height = axes_pos(4) * axes_pixel_ratio * image_ratio;
new_y = axes_pos(2) + (axes_pos(4) - new_height)/2;
new_pos = [ axes_pos(1) new_y axes_pos(3) new_height ];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Color Panel            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------
% --- Executes on slider movement.
function slider_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to slider_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gamma = get(hObject, 'Value');
axes(handles.axes_image);
cmap = imadjust(gray, [ 0 1 ], [], gamma);
colormap(cmap);
% Must reset the current axes to axes_overlay to get mouse clicks
axes(handles.axes_overlay);

%--------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function slider_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Main Panel            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------
% --- Executes on button press in but_done.
function but_done_Callback(hObject, eventdata, handles)
% hObject    handle to but_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
y_dim = handles.image_data.y_dim;
z_dim = handles.image_data.z_dim;

if ~isempty(handles.fig_data.regions)
    rect_bounds = zeros(size(handles.fig_data.rect_pos));
    rect_bounds(:,1) = round(handles.fig_data.rect_pos(:,1).*y_dim + 1);
    rect_bounds(:,2) = round((handles.fig_data.rect_pos(:,1) + ...
        handles.fig_data.rect_pos(:,3)).*y_dim);
    rect_bounds(:,3) = round((1 - (handles.fig_data.rect_pos(:,2) + ...
        handles.fig_data.rect_pos(:,4))).*z_dim + 1);
    rect_bounds(:,4) = round((1 - handles.fig_data.rect_pos(:,2)).*z_dim);
    handles.output = rect_bounds;
else
    handles.output = [];
end

guidata(hObject, handles);
uiresume;



%--------------------------------------------------------------------
% --- Executes on button press in but_cancel.
function but_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to but_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = [];
guidata(hObject, handles);
uiresume;

