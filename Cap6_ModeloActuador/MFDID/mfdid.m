function varargout = mfdid(varargin)
%MFDID M-file for mfdid.fig
%      MFDID, by itself, creates a new MFDID or raises the existing
%      singleton*.
%
%      H = MFDID returns the handle to a new MFDID or the handle to
%      the existing singleton*.
%
%      MFDID('Property','Value',...) creates a new MFDID using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to mfdid_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MFDID('CALLBACK') and MFDID('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MFDID.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mfdid

% Last Modified by GUIDE v2.5 15-Feb-2005 09:12:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mfdid_OpeningFcn, ...
                   'gui_OutputFcn',  @mfdid_OutputFcn, ...
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


% --- Executes just before mfdid is made visible.
function mfdid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

	handles.ind_ck = 0;
	if handles.ind_ck == 0,
		warning off
		echo off
	end

%	fprintf(1,'\n%d\n',nargin);
%	fprintf(1,'\n%s\n',varargin{1});

	if nargin > 3,
		if (varargin{1} == '3') | (varargin{1} == 'ABT1'),
			try,
%		        load ABT_exp02
		        load ABT_exp01
				ind_input = 3;
			catch,
				ind_input = 0;
			end
		elseif (varargin{1} == '2') | (varargin{1} == 'ABT0'),
			try,
		        load H_ABT
				ind_input = 2;
			catch,
				ind_input = 0;
			end
		elseif (varargin{1} == '1') | (varargin{1} == 'BIS'),
			try,
		        load Hyf_nm01
				ind_input = 1;
			catch,
				ind_input = 0;
			end
		else,
			ind_input = 0;
		end
	else,
		try,
	        load Hyf_nm01
			ind_input = 1;
		catch,
			try,
		        load H_ABT
				ind_input = 2;
			catch,
				try,
%			        load ABT_exp02
			        load ABT_exp01
					ind_input = 3;
				catch,
					ind_input = 0;
				end
			end
		end
	end

    try,
    	if ind_input == 1,

	        handles.freq = freq;					clear freq;
	        handles.H_exp = Hyf;					clear Hyf

	        handles.n_o = 2;
			handles.n_i = 2;

			handles.n_io = handles.n_i*handles.n_o;

	        handles.n_pole = 8;
	        handles.n_zero = 8*ones(handles.n_io,1);

	    	handles.ind_SR_o = [6,6].';
	    	handles.ind_SR_i = [2,1].';

	    	Sig_H = sigH;							clear sigH
	    	Sig_H_inv = sigH_inv;					clear sigH_inv

    	elseif ind_input == 2,

	        handles.freq = freq;					clear freq;
	        handles.H_exp = Hyf;					clear Hyf

	        handles.n_o = 2;
			handles.n_i = 2;

			handles.n_io = handles.n_i*handles.n_o;

			%
			%	ABT System
			%	ys1, ys2, yAMD
			%	yga, udr
			%

	        handles.n_pole = 7;
	        handles.n_zero = [7;7;5;6];

	    	handles.ind_SR_o = [6,6].';
	    	handles.ind_SR_i = [2,1].';

			coh = Cyf;								clear Cyf

    	elseif ind_input == 3,

	        handles.freq = freq;					clear freq;
	        handles.H_exp = Hyf;					clear Hyf

	        handles.n_o = 3;
			handles.n_i = 2;

			handles.n_io = handles.n_i*handles.n_o;

			%
			%	ABT System
			%	ys1, ys2, yAMD
			%	yga, udr
			%

	        handles.n_pole = 7;
	        handles.n_zero = [7;7;5;5;6;4];

	    	handles.ind_SR_o = [6,6,0].';
	    	handles.ind_SR_i = [2,1].';

			coh = Cyf;								clear Cyf

		else,

	        handles.freq = [1];
	        handles.H_exp = [1];

	        handles.n_o = 1;
			handles.n_i = 1;

			handles.n_io = handles.n_i*handles.n_o;

	        handles.n_pole = 1;
	        handles.n_zero = 1;

	    	handles.ind_SR_o = [6];
	    	handles.ind_SR_i = [0];

			coh = [1];

    	end
    catch,
        handles.freq = [1];
        handles.H_exp = [1];

        handles.n_o = 1;
		handles.n_i = 1;

		handles.n_io = handles.n_i*handles.n_o;

        handles.n_pole = 1;
        handles.n_zero = 1;

    	handles.ind_SR_o = [6];
    	handles.ind_SR_i = [0];

		coh = [1];
    end

	%	Weighting Function
        try,
	        handles.W_H = (1./Sig_H).^(1/2);	clear Sig_H
        catch,
			try,
		        handles.W_H = 1 ./ ...
		        	(sqrt(1-coh.^2) ./ abs(coh) ...
		        		.* abs(handles.H_exp));
		    catch,
		        handles.W_H = ones(size(handles.H_exp));
		    end
        end

		if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
		end

		if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
		end


        try,
	        handles.W_H_inv = (1./Sig_H_inv).^(1/2);
	        									clear Sig_H_inv
		catch,
			try,
		        handles.W_H_inv = 1 ./ ...
		        	(sqrt(1-coh.^2) ./ abs(coh) ...
		        		./ abs(handles.H_exp));
        	catch,
        		handles.W_H_inv = ones(size(handles.H_exp));
        	end
        end
		clear coh

		if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.W_H_inv).')
		[th1,th2] = view;
		v  = [sin(th1*pi/180),-cos(th1*pi/180),sin(th2*pi/180)];
		ax = [1,0,0];
		ax - dot(v / norm(v),ax) * v / norm(v);
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H_{inv} ', ...
			num2str(sum(sum(~isfinite(handles.W_H_inv))))])
		pause
		end

		if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.W_H_inv).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H_{inv} ', ...
			num2str(sum(sum(~isfinite(handles.W_H_inv))))])
		pause
		end

	if size(handles.H_exp,2) ~= handles.n_io,
		warning('Wrong Sensor Size!')
	end

	handles = data_pre(handles,0);

	% Choose default command line output for mfdid
	handles.output = hObject;

    % Update handles structure
	guidata(hObject, handles);

	% This sets up the initial plot - only do when we are invisible
	% so window can get raised using mfdid.
	if strcmp(get(hObject,'Visible'),'off')
	end

	% UIWAIT makes mfdid wait for user response (see UIRESUME)
	% uiwait(handles.fig_mfdid);


% --- Outputs from this function are returned to the command line.
function varargout = mfdid_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function fig_mfdid_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fig_mfdid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function handles = data_pre(handles,ind_con);

if ind_con == 0,
	%	Shape Check
		ind_sc = ~isfinite(handles.H_exp);
		handles.H_exp(ind_sc) = 0;
		handles.W_H(ind_sc) = 0;
		handles.W_H(~isfinite(handles.W_H)) = 0;

		ind_sc = ~isfinite(1./ handles.H_exp);
		handles.W_H_inv(ind_sc) = 0;
		handles.W_H_inv(~isfinite(handles.W_H_inv)) = 0;

    %   DC Terms
		ind_dc = handles.freq <= 0;
		handles.H_exp(ind_dc,:) = [];
		handles.W_H(ind_dc,:) = [];
		handles.W_H_inv(ind_dc,:) = [];
		handles.freq(ind_dc,:) = [];

    %   Scaling
        handles.f_eval	= max(handles.freq) / 3;
		set(handles.edit_f_eval,'String',...
	    	{num2str(handles.f_eval)});

    %   Frequency Range
        handles.ind_freq = ...
        	handles.freq >= min(handles.freq) & ...
        	handles.freq <= max(max(handles.freq)*0.65, ...
        		min(handles.freq));

	%
	%	Intial System Model
	%
	    handles.ind_sys = 0;
	    handles.sys_pool = {};

	%
	%	Fixed Pole Zero Model
	%
		handles.fixed_pole = [];
		for ii = 1:handles.n_io,
			handles.fixed_zero{ii} = [];
		end
end

	%
	%	Initial Screen
	%
	    handles.ind_f = 1;

	    axes(handles.axes_mag) % Select the proper axes
	    semilogy(handles.freq,abs(handles.H_exp));
	    xlabel('Frequency (Hz)')
	    ylabel('Magnitude')
		grid on
		set(gca,'YMinorGrid','off')

	    axes(handles.axes_pha) % Select the proper axes
	    plot(handles.freq,my_angle(handles.H_exp));
	    xlabel('Frequency (Hz)')
	    ylabel('Phase (^o)')
		grid on

		set(handles.edit_n_Poles,'String',...
	    	{num2str(handles.n_pole)});
		set(handles.edit_n_Zeros,'String',...
	    	{num2str(handles.n_zero(handles.ind_f))});


function [eV,Del_eV]=mrtf_lm(theta, ...
    H_exp, W_H, H_exp_inv, W_H_inv, ...
	n_io, n_pole_R, n_zero_R, ...
    ind_Real, ind_SR, ind_Inv, ...
    TBa, sE, ...
	Phi_a, Phi_a_c, Phi_B, Phi_B_c)
%MRTF_LM  MIMO Rational Transfer Function Estimator
%	for the structural system.
%
%	Levenberg-Marquardt Method
%	Stacked Matrix Form
%
%	Kim, Saang Bum; 8/18/2003 9:16PM
%	Copyright (c) 1992-2003 by SeoHa Lab.
%	$Revision: 0.0.0.1 $ $Date: 12/30/01 6:00PM $

%	Basic Environment Setting

	global eN_pool

	%	Parameters
	theta_a = [1; ...
		theta(1:n_pole_R)];
   	theta_a = apolystab(theta_a,ind_Real);
    theta_B = theta(n_pole_R+1:end);

	theta_a(1) = [];
	theta = [theta_a;theta_B];

	%	Error Norm Check
       if ind_SR == 0,
		H_est = (Phi_B * theta_B) ...
			./ (sE.^n_pole_R + Phi_a * theta_a);
		else,
			H_est = (Phi_B_c * theta_B ...
				+ Phi_a_c * theta_a) ...
				./ (sE.^n_pole_R + Phi_a * theta_a);
		end

		eV = (H_exp-H_est).*W_H;

		if ind_Inv == 1,
			H_est_inv = 1 ./ H_est;
			eV = [eV;
				(H_exp_inv-H_est_inv).*W_H_inv];
		end

		ind_ck = ~isfinite(eV);
%		ind_ck = logical(sign(sum(ind_ck.')));
		eV(ind_ck) = 0;

		eN = eV'*eV;
		eN_pool = [eN_pool;eN];

    %	Gradient
	if ind_SR == 0,
		ind_K1 = size(Phi_a,2);
		ind_K2 = size(Phi_B,2);

		Del1_eV = ((H_est ./ (sE.^n_pole_R + Phi_a * theta_a)) ...
			*ones(1,ind_K1)) .* Phi_a;

		Del2_eV = - Phi_B ./ ...
			((sE.^n_pole_R + Phi_a * theta_a)*ones(1,ind_K2));

		Del_eV = (W_H*ones(1,ind_K1 + ind_K2)) ...
			.*[Del1_eV Del2_eV];
	else,
		ind_K1 = size(Phi_a,2);
		ind_K2 = size(Phi_B_c,2);

		Del1_eV = ((H_est ./ (sE.^n_pole_R + Phi_a * theta_a)) ...
			*ones(1,ind_K1)) .* Phi_a ...
			- Phi_a_c ./ ...
            ((sE.^n_pole_R + Phi_a * theta_a)*ones(1,ind_K1));

			Del2_eV = - Phi_B_c ./ ...
				((sE.^n_pole_R + Phi_a * theta_a)*ones(1,ind_K2));

			Del_eV = (W_H*ones(1,ind_K1 + ind_K2)) ...
				.*[Del1_eV Del2_eV];
	end

	if ind_Inv == 1,
		if ind_SR == 0,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B,2);

			Del1_eV_inv = -Phi_a ...
				./ (Phi_B * theta_B *ones(1,ind_K1));

			Del2_eV_inv = ...
				(( ...
					H_est_inv ...
					./ ...
					(Phi_B * theta_B) ...
				)*ones(1,ind_K2)) ...
				.* Phi_B;

			Del_eV_inv = (W_H_inv*ones(1,ind_K1 + ind_K2)) ...
				.*[Del1_eV_inv Del2_eV_inv];
		else,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B_c,2);

			Del1_eV_inv = -Phi_a ...
				./ ((Phi_B_c * theta_B + Phi_a_c * theta_a) ...
					*ones(1,ind_K1)) ...
				+ ( ...
					(H_est_inv ...
					./ ...
					(Phi_B_c * theta_B + Phi_a_c * theta_a)) ...
					*ones(1,ind_K1)) ...
					.* Phi_a_c;

			Del2_eV_inv = ( ...
				(H_est_inv ...
				./ ...
				(Phi_B_c * theta_B +  Phi_a_c * theta_a)) ...
				* ones(1,ind_K2)) ...
				.* Phi_B_c;

			Del_eV_inv = (W_H_inv*ones(1,ind_K1 + ind_K2)) ...
				.*[Del1_eV_inv Del2_eV_inv];
		end
		Del_eV = [Del_eV;
			Del_eV_inv];
	end

	ind_ck = ~isfinite(Del_eV);
	ind_ck = logical(sign(sum(ind_ck.')));
	Del_eV(ind_ck,:) = 0;

	eV = [real(eV);imag(eV)];
	Del_eV = [real(Del_eV);imag(Del_eV)];

return

%
%	FINE
%


function [phase_H]=my_angle(Hyx)
%my_angle  Give a phase angle between -180 and 180 degree
%
%	[phase_H] = my_angle(Hyx)

%	Kim, Saang Bum; 2001. 4. 26. Thu. pm 5:2:53
%	Copyright (c) 1992-2003 by SeoHa Lab.
%	$Revision: 0.0.0.1 $ $Date: 7/14/2003 11:34AM $

	phase_H = mod(angle(Hyx)*180/pi+180,360)-180;

%
%	FINE
%


% --- Executes on button press in pushbutton_Prev_IO.
function pushbutton_Prev_IO_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Prev_IO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.ind_f = mod(handles.ind_f - 2,handles.n_io) + 1;

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in pushbutton_Next_IO.
function pushbutton_Next_IO_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Next_IO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.ind_f = mod(handles.ind_f,handles.n_io) + 1;

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in pushbutton_LLS.
function pushbutton_LLS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	handles.ind_ck = 0;
	if handles.ind_ck == 0,
		warning off
		echo off
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end


%
%	PART I.		Basic Parameter Setting
%
	im = sqrt(-1);

	freq     = handles.freq(handles.ind_freq);
	H_exp    = handles.H_exp(handles.ind_freq,:);
	H_exp_inv= 1./H_exp;
	W_H      = handles.W_H(handles.ind_freq,:);
%	W_H      = ones(size(handles.W_H(handles.ind_freq,:)));
	W_H_inv  = handles.W_H_inv(handles.ind_freq,:);
%	W_H_inv  = ones(size(handles.W_H_inv(handles.ind_freq,:)));
	n_o      = handles.n_o;
	n_i      = handles.n_i;
	n_io     = handles.n_io;
	ind_SR_o = handles.ind_SR_o;
	ind_SR_i = handles.ind_SR_i;
	n_pole	 = handles.n_pole;
	n_zero	 = handles.n_zero;
	n_pole_F = 0;
	n_zero_F = zeros(n_io,1);
    n_pole_R = n_pole - n_pole_F;
    n_zero_R = n_zero - n_zero_F;
	pF		 = handles.fixed_pole;
	zF		 = handles.fixed_zero;

	s = im*2*pi*freq;
    n_freq = length(freq);

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(H_exp))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(W_H))))])
		pause
	end

	%
	%	Option : Real Flag
	%
		realFlag = 0;
	    if get(handles.togglebutton_Real,'Value') > 0,
		    realFlag = 1;
		end

	%
	%	Option : Fixed
	%
		ind_Fixed = 0;
	    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
		    ind_Fixed = 1;

			for ii=1:length(pF),
				H_exp = H_exp ...
					.* ...
					((s - pF(ii)) ...
						*ones(1,n_io));
				W_H = W_H ...
					./ ...
					((s - pF(ii)) ...
						*ones(1,n_io));
			end
			for jj=1:n_io,
			for ii=1:length(zF{jj}),
				H_exp(:,jj) = H_exp(:,jj) ...
					./ (s - zF{jj}(ii));
				W_H(:,jj) = W_H(:,jj) ...
					.* (s - zF{jj}(ii));
			end
			end

		    ind_sc = ~isfinite(H_exp);
		    H_exp(ind_sc) = 0;
			W_H(ind_sc) = 0;

		    n_pole_F = length(pF);
			for ii=1:n_io,
				n_zero_F(ii) = length(zF{ii});
		    end
		    n_pole_R = n_pole - n_pole_F;
		    n_zero_R = n_zero - n_zero_F;
	    end

	%
	%	Option : Structural Relationships
	%
		ind_SR = 0;
	    TBa = cell(n_io,1);
	    if get(handles.togglebutton_SR,'Value') > 0,
		    ind_SR = 1;

		    for ii=1:n_i,
			    for jj=1:n_o,
			    	switch ind_SR_i(ii),
			    		case 0,
			    			TBa{jj+(ii-1)*n_o} = [];
			    		case 1,
			    			switch ind_SR_o(jj),
			    				case {0,1,4},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2,5},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
                                case {6},
                                    TBa{jj+(ii-1)*n_o} = zeros(4,4);
			    			end
			    		case 2,
			    			switch ind_SR_o(jj),
			    				case {0,1},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
			    				case {4,5,6},
			    					TBa{jj+(ii-1)*n_o} = eye(2,2);
			    			end
			    	end
				end
			end

			if ind_Fixed == 1,
				Phi_pole = prod(-pF);
				Sig_pole = 0;
				for ii=1:length(pF),
					fixed_pole = pF;
					fixed_pole(ii) = [];
					Sig_pole = Sig_pole + prod(-fixed_pole);
				end

			    for ii=1:n_i,
				    for jj=1:n_o,

						if ((ind_SR_i(ii) == 2) & (ind_SR_o(jj) > 3)),
					    	fixed_zero = zF{jj+(ii-1)*n_o};
							Phi_zero = ...
								prod(-fixed_zero);
							Sig_zero = 0;

							if Phi_zero == 0,
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
							else,
								for kk = 1:length(fixed_zero),
									fixed_zero_kk  = fixed_zero;
									fixed_zero_kk(kk) = [];
									Sig_zero = Sig_zero ...
										+ prod(-fixed_zero_kk);
								end
		    					TBa{jj+(ii-1)*n_o} = ...
									[Phi_pole/Phi_zero, ...
										( ...
											Sig_pole ...
											-Phi_pole/Phi_zero*Sig_zero ...
										) / Phi_zero;
									0,Phi_pole/Phi_zero];
							end
						end

					end
				end
			end
	    end

	%
	%	Option : Inverse Model
	%
	    ind_Inv = 0;
	    if get(handles.togglebutton_Inv,'Value') > 0,
		    ind_Inv = 1;

			if ind_Fixed == 1,
				for ii=1:length(pF),
					H_exp_inv = H_exp_inv ...
						./ ...
						( ...
							(s - pF(ii)) ...
							*ones(1,n_io) ...
						);

					W_H_inv = W_H_inv ...
						.* ...
						( ...
							(s - pF(ii)) ...
							*ones(1,n_io) ...
						);
				end
				for jj=1:n_io,
					for ii=1:length(zF{jj}),
						W_H_inv(:,jj) = W_H_inv(:,jj) ...
							./ (s - zF{jj}(ii));
						H_exp_inv(:,jj) = H_exp_inv(:,jj) ...
							.* (s - zF{jj}(ii));
					end
				end
			end

			ind_sc = ~isfinite(H_exp_inv);
			H_exp_inv(ind_sc) = 0;
			W_H_inv(ind_sc) = 0;
	    end


%
%	PART II.	Pre Process
%				Transform to Stacked Matrices
%
	H_exp = reshape(H_exp,n_freq*n_io,1);
	W_H	  = reshape(W_H  ,n_freq*n_io,1);

	if ind_Inv == 1,
		H_exp_inv = reshape(H_exp_inv,n_freq*n_io,1);
		W_H_inv	  = reshape(W_H_inv  ,n_freq*n_io,1);
	end

	sE = [];
	for ii=1:n_io,
		sE = [sE; ...
			s];
	end

	n_max = max(n_pole_R, max(n_zero_R) + 1);
	Phi = ones(n_freq,n_max);
	for ii = 2:n_max,
		Phi(:,ii) = s.^(ii-1);
	end

	%
	%	ARMA Parameters
	%
		%
		%	Phi_a
		%
			Phi_a = zeros(n_freq*n_io,n_pole_R);
			for ii = 1:n_io,
				Phi_a([1:n_freq]+(ii-1)*n_freq,:) = ...
					Phi(:,[n_pole_R:-1:1]);
			end

		%
		%	Phi_B
		%
			Phi_B = zeros(n_freq*n_io,sum(n_zero_R)+n_io);
			for ii=1:n_io,
				ind_B = (n_zero_R(ii)+1):-1:1;
				Phi_B([1:n_freq]+(ii-1)*n_freq, ...
					[1:(n_zero_R(ii)+1)] + ...
					ii-1+sum(n_zero_R(1:ii-1))) = ...
					Phi(:,ind_B);
			end

		if ind_SR == 1,
		%
		%	Phi_a_c & Phi_B_c
		%
			Phi_a_c = zeros(n_freq*n_io,n_pole_R);
			Phi_B_c = Phi_B;
			for ii=n_io:-1:1,
				n_c = size(TBa{ii},2);
				Phi_a_c([1:n_freq]+(ii-1)*n_freq, ...
					(n_pole_R+1-[n_c:-1:1])) = ...
					Phi(:,[n_c:-1:1]) * TBa{ii};
				Phi_B_c(:,n_zero_R(ii)+2+ ...
					sum(n_zero_R(1:ii-1))+(ii-1)- ...
					[n_c:-1:1]) = [];
			end
		end


%
%	PART III.	Main Process
%				Linear Least Squares Estimation
%
	%
	%	Weighting Adjustment
	%
%		W_H = ones(size(W_H));
		W_H = W_H;
%		W_H = abs(H_exp);
%		W_H = W_H .* abs(H_exp);

	%
	%	Force & Stiffness
	%
		f = W_H .* H_exp .* (sE.^n_pole_R);
%		f = H_exp .* (sE.^n_pole_R) .* W_H;

		if ind_SR == 0,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B,2);
			K = (W_H * ones(1, ind_K1+ind_K2)) .* ...
				[-(H_exp*ones(1,ind_K1)) .* Phi_a, ...
					Phi_B];
%			K = [-(H_exp*ones(1,ind_K1)) .* Phi_a, ...
%					Phi_B] ...
%				.* (W_H * ones(1, ind_K1+ind_K2));
		else,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B_c,2);
			K = (W_H * ones(1, ind_K1+ind_K2)) .* ...
				[-(H_exp*ones(1,ind_K1)) .* Phi_a + Phi_a_c, ...
					Phi_B_c];
%			K = [-(H_exp*ones(1,ind_K1)) .* Phi_a + Phi_a_c, ...
%					Phi_B_c] ...
%				.* (W_H * ones(1, ind_K1+ind_K2));
		end
	if ind_Inv == 1,
		f = [f;
			W_H_inv .* (sE.^n_pole_R)];

		if ind_SR == 0,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B,2);
			K = [K;
					(W_H_inv * ones(1, ind_K1+ind_K2)) .* ...
					[-Phi_a, ...
						(H_exp_inv*ones(1,ind_K2)) .* Phi_B];
				];
		else,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B_c,2);
			K = [K;
					(W_H_inv * ones(1, ind_K1+ind_K2)) .* ...
					[-Phi_a + ...
						(H_exp_inv*ones(1,ind_K1)) .* Phi_a_c, ...
						(H_exp_inv*ones(1,ind_K2)) .* Phi_B_c];
				];
		end
	end

	%
	%	Displacement
	%
	if realFlag == 0,
		theta = (K'*K) \ (K'*f);

		if handles.ind_ck > 0,
			figure(handles.ind_ck)
			clf
			bar(abs(real(theta) ./ abs(theta)))
			hold on
			bar(-abs(imag(theta) ./ abs(theta)),'r')
			hold off
		end

	else,
%		theta = (real(K'*K) \ real(K'*f) + ...
%			imag(K'*K) \ imag(K'*f))/2;
		theta = real(K'*K) \ real(K'*f);
	end


%
%	PART IV.	Post Process
%
	%
	%	Parameter Vectors
	%
		theta_a_R = [1;theta(1:n_pole_R)];
    	theta_a_R = apolystab(theta_a_R,realFlag);
    	theta_B_R = cell(n_io,1);
		if ind_SR == 0,
			for jj=1:n_i,
				for ii=1:n_o,
					theta_B_R{ii+(jj-1)*n_o} = ...
						theta([1:(n_zero_R(ii+(jj-1)*n_o)+1)] + ...
							sum(n_zero_R(1:(ii-1+(jj-1)*n_o))) + ...
							(ii-1+(jj-1)*n_o) + ...
							n_pole_R);

				end
			end
		else,
			ind_th = n_pole_R;
			for jj=1:n_i,
				for ii=1:n_o,
					n_c = size(TBa{ii+(jj-1)*n_o},2);
					theta_B_R{ii+(jj-1)*n_o} = ...
						[theta([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c] + ...
							ind_th);
							TBa{ii+(jj-1)*n_o}*...
								theta_a_R(n_pole_R+2-[n_c:-1:1]);
						];
					ind_th = ind_th + ...
						n_zero_R(ii+(jj-1)*n_o)+1-n_c;
				end
			end
		end

	%
	%	Rational Polynomial Model
	%
		for jj=1:n_i,
			for ii=1:n_o,
				num{ii,jj} = theta_B_R{ii+(jj-1)*n_o}.';
				den{ii,jj} = theta_a_R.';
			end
		end
		sys_R = tf(num,den);

		if ind_Fixed == 1,
			[z,p,k] = zpkdata(sys_R);
			for jj=1:n_i,
				for ii=1:n_o,
	        		p{ii,jj} = [p{ii,jj}; ...
	                    handles.fixed_pole];
					z{ii,jj} = [z{ii,jj}; ...
	                        handles.fixed_zero{ii+(jj-1)*n_o}];
				end
			end
			sys_R = zpk(z,p,k);
		else,
			handles.fixed_pole = [];
			for ii = 1:handles.n_io,
				handles.fixed_zero{ii} = [];
			end
		end

	%
	%	Data Base & Graphic
	%
	    handles.ind_sys = length(handles.sys_pool) + 1;
	    handles.sys_pool{handles.ind_sys} = ...
	 		sys_R;

		str_ind_sys = [ ...
			num2str(handles.ind_sys), ...
			'/', ...
			num2str(length(handles.sys_pool))];
		set(handles.text_ind_sys, 'String', str_ind_sys);
		n_error_norm(handles);

		draw_n(handles);
		draw_frf(handles);
	    try,
		    pz_window(handles);
		end

    %   Update handles structure
		guidata(hObject, handles);


function a = apolystab(a,realFlag)
%APOLYSTAB  Stabilize filter, analog
%   inputs: a - denominator polynomial
%           realFlag - 1 for real, 0 for complex
%   returns stabilized denoninator polynomial
if length(a)>0
    v=roots(a);
    vind=find(real(v)>0);
    v(vind)=-v(vind);
    a=poly(v);
    if realFlag
        a=real(a);
    end
    a = a.';
end


function n_error_norm(handles)

	H_est = freqresp(handles.sys_pool{handles.ind_sys}, ...
		handles.freq(handles.ind_freq)*2*pi);
	H_mdl = [];
	for ii=1:handles.n_i,
		for jj=1:handles.n_o,
			H_mdl = [H_mdl, (squeeze(H_est(jj,ii,:)))];
		end
	end

	H_exp = handles.H_exp(handles.ind_freq,:);
	W_H = handles.W_H(handles.ind_freq,:);

	eB = H_exp.*W_H;
	eV = (H_exp - H_mdl).*W_H;

    [m,n] = size(eB);
	eB = reshape(eB,m*n,1);
	eV = reshape(eV,m*n,1);

    eB = eB'*eB;
    eV = eV'*eV;

	set(handles.text_eN,'String',num2str(eV));
	set(handles.text_eB,'String',num2str(eB));


function edit_n_Poles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_Poles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_Poles as text
%        str2double(get(hObject,'String')) returns contents of edit_n_Poles as a double

handles.n_pole = str2double(get(hObject,'string'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_n_Poles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_Poles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_n_Zeros_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_Zeros (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_Zeros as text
%        str2double(get(hObject,'String')) returns contents of edit_n_Zeros as a double

handles.n_zero(handles.ind_f) = str2double(get(hObject,'string'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_n_Zeros_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_Zeros (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_Pole_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Pole_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Pole_R as text
%        str2double(get(hObject,'String')) returns contents of edit_Pole_R as a double


% --- Executes during object creation, after setting all properties.
function edit_Pole_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Pole_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_Pole_I_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Pole_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Pole_I as text
%        str2double(get(hObject,'String')) returns contents of edit_Pole_I as a double


% --- Executes during object creation, after setting all properties.
function edit_Pole_I_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Pole_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_Zero_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Zero_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zero_R as text
%        str2double(get(hObject,'String')) returns contents of edit_Zero_R as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
%    errordlg('You must enter a numeric value','Bad Input','modal')
end
% proceed with callback...


% --- Executes during object creation, after setting all properties.
function edit_Zero_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Zero_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_Zero_I_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Zero_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zero_I as text
%        str2double(get(hObject,'String')) returns contents of edit_Zero_I as a double


% --- Executes during object creation, after setting all properties.
function edit_Zero_I_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Zero_I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton_Close.
function pushbutton_Close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.fig_mfdid)


% --- Executes on button press in pushbutton_SM.
function pushbutton_SM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	handles.ind_ck = 0;
	if handles.ind_ck == 0,
		warning off
		echo off
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end


%
%	PART I.		Basic Parameter Setting
%
	im = sqrt(-1);

	freq     = handles.freq(handles.ind_freq);
	H_exp    = handles.H_exp(handles.ind_freq,:);
	H_exp_inv= 1./H_exp;
	n_o      = handles.n_o;
	n_i      = handles.n_i;
	n_io     = handles.n_io;
	ind_SR_o = handles.ind_SR_o;
	ind_SR_i = handles.ind_SR_i;
	W_H      = handles.W_H(handles.ind_freq,:);
%	W_H      = ones(size(handles.W_H(handles.ind_freq,:)));
	W_H_inv  = handles.W_H_inv(handles.ind_freq,:);
%	W_H_inv  = ones(size(handles.W_H_inv(handles.ind_freq,:)));
	n_pole	 = handles.n_pole;
	n_zero	 = handles.n_zero;
	n_pole_F = 0;
	n_zero_F = zeros(n_io,1);
    n_pole_R = n_pole - n_pole_F;
    n_zero_R = n_zero - n_zero_F;
	pF		 = handles.fixed_pole;
	zF		 = handles.fixed_zero;

	s = im*2*pi*freq;
    n_freq = length(freq);

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(H_exp))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(W_H))))])
		pause
	end

	%
	%	Option : Real Flag
	%
		realFlag = 0;
	    if get(handles.togglebutton_Real,'Value') > 0,
		    realFlag = 1;
		end

	%
	%	Option : Fixed
	%
		ind_Fixed = 0;
	    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
		    ind_Fixed = 1;

			for ii=1:length(pF),
				H_exp = H_exp ...
					.* ...
					((s - pF(ii)) ...
						*ones(1,n_io));
				W_H = W_H ...
					./ ...
					((s - pF(ii)) ...
						*ones(1,n_io));
			end
			for jj=1:n_io,
			for ii=1:length(zF{jj}),
				H_exp(:,jj) = H_exp(:,jj) ...
					./ (s - zF{jj}(ii));
				W_H(:,jj) = W_H(:,jj) ...
					.* (s - zF{jj}(ii));
			end
			end

		    ind_sc = ~isfinite(H_exp);
		    H_exp(ind_sc) = 0;
			W_H(ind_sc) = 0;

		    n_pole_F = length(pF);
			for ii=1:n_io,
				n_zero_F(ii) = length(zF{ii});
		    end
		    n_pole_R = n_pole - n_pole_F;
		    n_zero_R = n_zero - n_zero_F;
	    end

	%
	%	Option : Structural Relationships
	%
		ind_SR = 0;
	    TBa = cell(n_io,1);
	    if get(handles.togglebutton_SR,'Value') > 0,
		    ind_SR = 1;

		    for ii=1:n_i,
			    for jj=1:n_o,
			    	switch ind_SR_i(ii),
			    		case 0,
			    			TBa{jj+(ii-1)*n_o} = [];
			    		case 1,
			    			switch ind_SR_o(jj),
			    				case {0,1,4},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2,5},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
                                case {6},
                                    TBa{jj+(ii-1)*n_o} = zeros(4,4);
			    			end
			    		case 2,
			    			switch ind_SR_o(jj),
			    				case {0,1},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
			    				case {4,5,6},
			    					TBa{jj+(ii-1)*n_o} = eye(2,2);
			    			end
			    	end
				end
			end

			if ind_Fixed == 1,
				Phi_pole = prod(-pF);
				Sig_pole = 0;
				for ii=1:length(pF),
					fixed_pole = pF;
					fixed_pole(ii) = [];
					Sig_pole = Sig_pole + prod(-fixed_pole);
				end

			    for ii=1:n_i,
				    for jj=1:n_o,

						if ((ind_SR_i(ii) == 2) & (ind_SR_o(jj) > 3)),
					    	fixed_zero = zF{jj+(ii-1)*n_o};
							Phi_zero = ...
								prod(-fixed_zero);
							Sig_zero = 0;

							if Phi_zero == 0,
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
							else,
								for kk = 1:length(fixed_zero),
									fixed_zero_kk  = fixed_zero;
									fixed_zero_kk(kk) = [];
									Sig_zero = Sig_zero ...
										+ prod(-fixed_zero_kk);
								end
		    					TBa{jj+(ii-1)*n_o} = ...
									[Phi_pole/Phi_zero, ...
										( ...
											Sig_pole ...
											-Phi_pole/Phi_zero*Sig_zero ...
										) / Phi_zero;
									0,Phi_pole/Phi_zero];
							end
						end

					end
				end
			end
	    end

	%
	%	Option : Inverse Model
	%
	    ind_Inv = 0;
	    if get(handles.togglebutton_Inv,'Value') > 0,
		    ind_Inv = 1;

			if ind_Fixed == 1,
				for ii=1:length(pF),
					H_exp_inv = H_exp_inv ...
						./ ...
						( ...
							(s - pF(ii)) ...
							*ones(1,n_io) ...
						);

					W_H_inv = W_H_inv ...
						.* ...
						( ...
							(s - pF(ii)) ...
							*ones(1,n_io) ...
						);
				end
				for jj=1:n_io,
					for ii=1:length(zF{jj}),
						W_H_inv(:,jj) = W_H_inv(:,jj) ...
							./ (s - zF{jj}(ii));
						H_exp_inv(:,jj) = H_exp_inv(:,jj) ...
							.* (s - zF{jj}(ii));
					end
				end
			end

			ind_sc = ~isfinite(H_exp_inv);
			H_exp_inv(ind_sc) = 0;
			W_H_inv(ind_sc) = 0;
	    end


%
%	PART II.	Pre Process
%				Transform to Stacked Matrices
%
	H_exp = reshape(H_exp,n_freq*n_io,1);
	W_H	  = reshape(W_H  ,n_freq*n_io,1);

	if ind_Inv == 1,
		H_exp_inv = reshape(H_exp_inv,n_freq*n_io,1);
		W_H_inv	  = reshape(W_H_inv  ,n_freq*n_io,1);
	end

	sE = [];
	for ii=1:n_io,
		sE = [sE; ...
			s];
	end

	n_max = max(n_pole_R, max(n_zero_R) + 1);
	Phi = ones(n_freq,n_max);
	for ii = 2:n_max,
		Phi(:,ii) = s.^(ii-1);
	end

	%
	%	ARMA Parameters
	%
		%
		%	Phi_a
		%
			Phi_a = zeros(n_freq*n_io,n_pole_R);
			for ii = 1:n_io,
				Phi_a([1:n_freq]+(ii-1)*n_freq,:) = ...
					Phi(:,[n_pole_R:-1:1]);
			end

		%
		%	Phi_B
		%
			Phi_B = zeros(n_freq*n_io,sum(n_zero_R)+n_io);
			for ii=1:n_io,
				ind_B = (n_zero_R(ii)+1):-1:1;
				Phi_B([1:n_freq]+(ii-1)*n_freq, ...
					[1:(n_zero_R(ii)+1)] + ...
					ii-1+sum(n_zero_R(1:ii-1))) = ...
					Phi(:,ind_B);
			end

		if ind_SR == 1,
		%
		%	Phi_a_c & Phi_B_c
		%
			Phi_a_c = zeros(n_freq*n_io,n_pole_R);
			Phi_B_c = Phi_B;
			for ii=n_io:-1:1,
				n_c = size(TBa{ii},2);
				Phi_a_c([1:n_freq]+(ii-1)*n_freq, ...
					(n_pole_R+1-[n_c:-1:1])) = ...
					Phi(:,[n_c:-1:1]) * TBa{ii};
				Phi_B_c(:,n_zero_R(ii)+2+ ...
					sum(n_zero_R(1:ii-1))+(ii-1)- ...
					[n_c:-1:1]) = [];
			end
		end


%
%	PART III.	Main Process
%				Steiglitz-McBride
%
	%
	%	Prepare the Iteration
	%
		sys_R = handles.sys_pool{handles.ind_sys};
		sys_F = sys_R;

		if ind_Fixed == 1,
			for jj=1:handles.n_i,
				for ii=1:handles.n_o,
					sys_F(ii,jj) = zpk(zF{ii+(jj-1)*n_o},pF, ...
						1);
					sys_R(ii,jj) = minreal(sys_R(ii,jj) ...
						/ sys_F(ii,jj));
	            end
			end
		end
		[num,den] = tfdata(sys_R);

		theta_a = den{1,1}.';
	    theta_B = [];
		for jj=1:n_i,
			for ii=1:n_o,
	            ind_srt = length(theta_a) - n_pole_R;

	            if ind_srt > 1,
	            	if handles.ind_ck == 1,
	                disp('check coef of pole')
	                disp(theta_a(1:ind_srt-1));
	                end
	                theta_a = theta_a(ind_srt:end);
	            end

	            theta_B_imsi = num{ii,jj}.' / theta_a(1);

	            ind_srt = length(theta_B_imsi) ...
	                - n_zero_R(ii+(jj-1)*n_o);

	            if ind_srt > 1,
	            	if handles.ind_ck == 1,
	                disp('check coef of zero')
	                disp(theta_B_imsi(1:ind_srt-1));
	                end
	                theta_B_imsi = theta_B_imsi(ind_srt:end);
	            end

	            theta_B = [theta_B; ...
					theta_B_imsi];
			end
		end
		theta_a = theta_a / theta_a(1);

		%	Initial Parameters

			theta_a(1) = [];
			theta = [theta_a;theta_B];

		    W_H_n = W_H;
		    if ind_Inv == 1,
		        W_H_inv_n = W_H_inv;
		    end

		%	The initial estimate
			if n_pole_R == 0,
				H_est = (Phi_B * theta_B);
			else,
				H_est = (Phi_B * theta_B) ...
					./ (sE.^n_pole_R + Phi_a * theta_a);
			end

			eV = (H_exp-H_est).*W_H;

			if ind_Inv == 1,
				H_est_inv = 1 ./ H_est;
				eV = [eV;
					(H_exp_inv-H_est_inv).*W_H_inv];
			end

			ind_ck = ~isfinite(eV);
			ind_ck = logical(sign(sum(ind_ck.')));
			eV(ind_ck,:) = 0;

			eN_pool = [eV'*eV];

			h_f_sm = figure;
			semilogy(0,eN_pool(1),'or')
			xlabel('Iteration #')
			ylabel('Error Norm')
			title('Steglitz-McBride Method')
			hold on

		%	Convergence Conditions
			maxiter = 30;
		    tol = 0.01;

			tol_r = tol *ones(size(theta));
			tol_a = tol *ones(size(theta));
			tol_b = zeros(size(theta));
			tol_n = maxiter;

	%
	%	Iteration
	%
	for ind_iter=1:tol_n,

    	theta_O = theta;

		%	Weighting Update
		if n_pole_R ~= 0,
	    	W_H_n = W_H ./ (sE.^n_pole_R + Phi_a*theta_a);
			W_H_n(~isfinite(W_H_n)) = 0;
		end

		if ind_Inv == 1,
			if n_zero_R ~= 0,
				W_H_inv_n = W_H_inv ./ (Phi_B*theta_B);
				W_H_inv_n(~isfinite(W_H_inv_n)) = 0;
			end
		end

		%
		%	Force & Stiffness
		%
			f = W_H_n .* H_exp .* (sE.^n_pole_R);

			if ind_SR == 0,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B,2);
				K = (W_H_n * ones(1, ind_K1+ind_K2)) .* ...
					[-(H_exp*ones(1,ind_K1)) .* Phi_a, ...
						Phi_B];
			else,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B_c,2);
				K = (W_H_n * ones(1, ind_K1+ind_K2)) .* ...
					[-(H_exp*ones(1,ind_K1)) .* Phi_a + Phi_a_c, ...
						Phi_B_c];
			end
		if ind_Inv == 1,
			f = [f;
				W_H_inv_n .* (sE.^n_pole_R)];

			if ind_SR == 0,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B,2);
				K = [K;
						(W_H_inv_n * ones(1, ind_K1+ind_K2)) .* ...
						[-Phi_a, ...
							(H_exp_inv*ones(1,ind_K2)) .* Phi_B];
					];
			else,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B_c,2);
				K = [K;
						(W_H_inv_n * ones(1, ind_K1+ind_K2)) .* ...
						[-Phi_a + ...
							(H_exp_inv*ones(1,ind_K1)) .* Phi_a_c, ...
							(H_exp_inv*ones(1,ind_K2)) .* Phi_B_c];
					];
			end
		end

		%
		%	Displacement
		%
		if realFlag == 0,
			theta = (K'*K) \ (K'*f);

			if handles.ind_ck > 0,
				figure(handles.ind_ck)
				clf
				bar(abs(real(theta) ./ abs(theta)))
				hold on
				bar(-abs(imag(theta) ./ abs(theta)),'r')
				hold off
			end

		else,
	%		theta = (real(K'*K) \ real(K'*f) + ...
	%			imag(K'*K) \ imag(K'*f))/2;
			theta = real(K'*K) \ real(K'*f);
		end

		%
		%	Parameter Vectors
		%
			theta_a = [1;theta(1:n_pole_R)];
	    	theta_a = apolystab(theta_a,realFlag);

		    theta_B = [];
			if ind_SR == 0,
				for jj=1:n_i,
					for ii=1:n_o,
						theta_B = [theta_B; ...
							theta([1:(n_zero_R(ii+(jj-1)*n_o)+1)] + ...
								sum(n_zero_R(1:(ii-1+(jj-1)*n_o))) + ...
								(ii-1+(jj-1)*n_o) + ...
								n_pole_R)];
					end
				end
			else,
				ind_th = n_pole_R;
				for jj=1:n_i,
					for ii=1:n_o,
						n_c = size(TBa{ii+(jj-1)*n_o},2);
						theta_B = [theta_B; ...
							[theta([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c] + ...
								ind_th);
								TBa{ii+(jj-1)*n_o}*...
									theta_a(n_pole_R+2-[n_c:-1:1]);
							]];
						ind_th = ind_th + ...
							n_zero_R(ii+(jj-1)*n_o)+1-n_c;
					end
				end
			end

		%
			theta_a(1) = [];
			theta = [theta_a;theta_B];
		%

		%
		%	Error Norm Check
		%
			if n_pole_R == 0,
				H_est = (Phi_B * theta_B);
			else,
				H_est = (Phi_B * theta_B) ...
					./ (sE.^n_pole_R + Phi_a * theta_a);
			end

			eV = (H_exp-H_est).*W_H;

			if ind_Inv == 1,
				H_est_inv = 1./H_est;

				eV = [eV;
					(H_exp_inv-H_est_inv).*W_H_inv];
			end

			ind_ck = ~isfinite(eV);
			ind_ck = logical(sign(sum(ind_ck.')));
			eV(ind_ck,:) = 0;

			eN_pool = [eN_pool;eV'*eV];

			figure(h_f_sm);
			semilogy(ind_iter,eN_pool(ind_iter+1),'o')

    	%
		%	3. Tolerance Check
		%
			tol_b = max(tol_b,theta);
			if (abs((theta - theta_O)./theta_O) < tol_r) ...
				& (abs(theta - theta_O) < tol_a.*tol_b),
				break
			end
	end

	figure(h_f_sm);
	semilogy(ind_iter,eN_pool(end),'or')
	hold off


%
%	PART IV.	Post Process
%
	%
	%	Parameter Vectors
	%
		theta_a_R = [1;theta(1:n_pole_R)];

	    theta_B_R = cell(n_io,1);
		for jj=1:n_i,
			for ii=1:n_o,
				theta_B_R{ii+(jj-1)*n_o} = ...
					theta([1:(n_zero_R(ii+(jj-1)*n_o)+1)] + ...
						sum(n_zero_R(1:(ii-1+(jj-1)*n_o))) + ...
						(ii-1+(jj-1)*n_o) + ...
						n_pole_R);
			end
		end

	%
	%	Rational Polynomial Model
	%
		for jj=1:n_i,
			for ii=1:n_o,
				num{ii,jj} = theta_B_R{ii+(jj-1)*n_o}.';
				den{ii,jj} = theta_a_R.';
			end
		end
		sys_R = tf(num,den);

		if ind_Fixed == 1,
			[z,p,k] = zpkdata(sys_R);
			for jj=1:n_i,
				for ii=1:n_o,
	        		p{ii,jj} = [p{ii,jj}; ...
	                    handles.fixed_pole];
					z{ii,jj} = [z{ii,jj}; ...
	                        handles.fixed_zero{ii+(jj-1)*n_o}];
				end
			end
			sys_R = zpk(z,p,k);
		else,
			handles.fixed_pole = [];
			for ii = 1:handles.n_io,
				handles.fixed_zero{ii} = [];
			end
		end

	%
	%	Data Base & Graphic
	%
	    handles.ind_sys = length(handles.sys_pool) + 1;
	    handles.sys_pool{handles.ind_sys} = ...
	 		sys_R;

		str_ind_sys = [ ...
			num2str(handles.ind_sys), ...
			'/', ...
			num2str(length(handles.sys_pool))];
		set(handles.text_ind_sys, 'String', str_ind_sys);
		n_error_norm(handles);

		draw_n(handles);
		draw_frf(handles);
	    try,
		    pz_window(handles);
		end

    %   Update handles structure
		guidata(hObject, handles);


% --- Executes on button press in pushbutton_Add_Sys.
function pushbutton_Add_Sys_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Add_Sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	n_o  = handles.n_o;
	n_i  = handles.n_i;
	n_io = handles.n_io;

	s_p_r = cellstr(get(handles.edit_Pole_R,'string'));
	s_p_i = cellstr(get(handles.edit_Pole_I,'string'));
	s_z_r = cellstr(get(handles.edit_Zero_R,'string'));
	s_z_i = cellstr(get(handles.edit_Zero_I,'string'));

	ind_fixed_pole = zeros(length(s_p_r),1);
	for ii=1:length(s_p_r),
		if lower(s_p_r{ii}(1)) == 'f',
			s_p_r{ii}(1) = [];
			ind_fixed_pole(ii) = 1;
		end
		if lower(s_p_i{ii}(1)) == 'f',
			s_p_i{ii}(1) = [];
			ind_fixed_pole(ii) = 1;
		end
	end

	ind_fixed_zero = zeros(length(s_z_r),1);
	for ii=1:length(s_z_r),
		if lower(s_z_r{ii}(1)) == 'f',
			s_z_r{ii}(1) = [];
			ind_fixed_zero(ii) = 1;
		end
		if lower(s_z_i{ii}(1)) == 'f',
			s_z_i{ii}(1) = [];
			ind_fixed_zero(ii) = 1;
		end
	end

    if ~isempty(s_p_r)
    	p_r = str2double(s_p_r);
    	p_i = str2double(s_p_i);
    else
        p_r = [];
        p_i = [];
    end
    if ~isempty(s_z_r)
    	z_r = str2double(s_z_r);
    	z_i = str2double(s_z_i);
    else
        z_r = [];
        z_i = [];
    end

	ind = isfinite(p_r) & isfinite(p_i);
    if ~isempty(p_r),
        ind_fixed_pole = logical(ind_fixed_pole .* ind);
    end

    p = p_r(ind) + sqrt(-1)*p_i(ind);
	handles.fixed_pole = p_r(ind_fixed_pole) + ...
		sqrt(-1)*p_i(ind_fixed_pole);

	ind = isfinite(z_r) & isfinite(z_i);
    if ~isempty(z_r),
    	ind_fixed_zero = logical(ind_fixed_zero .* ind);
    end
	z = z_r(ind) + sqrt(-1)*z_i(ind);
	handles.fixed_zero{handles.ind_f} = z_r(ind_fixed_zero) + ...
		sqrt(-1)*z_i(ind_fixed_zero);

	handles.n_pole = length(p);
	handles.n_zero(handles.ind_f) = length(z);

	set(handles.edit_n_Poles,'String',...
    	{num2str(handles.n_pole)});
	set(handles.edit_n_Zeros,'String',...
    	{num2str(handles.n_zero(handles.ind_f))});

	[zN,pN,kN] = zpkdata(handles.sys_pool{handles.ind_sys});

	for jj=1:handles.n_i,
		for ii=1:handles.n_o,
			pN{ii,jj} = p;
		end
	end
	zN{handles.ind_f} = z;

    %   stabilizing

	%   Scaling
    ind_f_eval = sum(handles.freq <= handles.f_eval);
    set(handles.edit_f_eval, 'String', num2str(handles.f_eval));

    for ii=1:handles.n_io,
	    sys_est= zpk(zN{ii},pN{ii},1);
	    H_est  = freqresp(sys_est,handles.freq(ind_f_eval)*2*pi);
	    H_exp  = handles.H_exp(ind_f_eval,ii);

%	    kN(ii) = abs(h_ref)/abs(h_imsi);
	    kN(ii) = H_exp/H_est;
    end

	sys_est = zpk(zN,pN,kN);
	sys_R = sys_est;

	%
	%	Option : Structural Relationships
	%
		ind_SR_o = handles.ind_SR_o;
		ind_SR_i = handles.ind_SR_i;
		pF = handles.fixed_pole;
		zF = handles.fixed_zero;

	    if get(handles.togglebutton_SR,'Value') > 0,
			for jj=1:handles.n_i,
				for ii=1:handles.n_o,
					sys_F = zpk(zF{ii+(jj-1)*n_o},pF, ...
						1);
					sys_R(ii,jj) = minreal(sys_est(ii,jj) / sys_F);
                end
			end

		    TBa = cell(n_io,1);
			for ii=1:n_i,
			    for jj=1:n_o,
			    	switch ind_SR_i(ii),
			    		case 0,
			    			TBa{jj+(ii-1)*n_o} = [];
			    		case 1,
			    			switch ind_SR_o(jj),
			    				case {0,1,4},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2,5},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
                                case {6},
                                    TBa{jj+(ii-1)*n_o} = zeros(4,4);
			    			end
			    		case 2,
			    			switch ind_SR_o(jj),
			    				case {0,1},
			    					TBa{jj+(ii-1)*n_o} = [];
			    				case {2},
			    					TBa{jj+(ii-1)*n_o} = [0];
			    				case {3},
			    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
			    				case {4,5,6},
			    					TBa{jj+(ii-1)*n_o} = eye(2,2);
			    			end
			    	end
				end
			end

		    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
				Phi_pole = prod(-pF);
				Sig_pole = 0;
				for ii=1:length(pF),
					fixed_pole = pF;
					fixed_pole(ii) = [];
					Sig_pole = Sig_pole + prod(-fixed_pole);
				end

			    for ii=1:n_i,
				    for jj=1:n_o,

						if ((ind_SR_i(ii) == 2) & (ind_SR_o(jj) > 3)),
					    	fixed_zero = zF{jj+(ii-1)*n_o};
							Phi_zero = ...
								prod(-fixed_zero);
							Sig_zero = 0;

							if Phi_zero == 0,
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
							else,
								for kk = 1:length(fixed_zero),
									fixed_zero_kk  = fixed_zero;
									fixed_zero_kk(kk) = [];
									Sig_zero = Sig_zero ...
										+ prod(-fixed_zero_kk);
								end
		    					TBa{jj+(ii-1)*n_o} = ...
									[Phi_pole/Phi_zero, ...
										( ...
											Sig_pole ...
											-Phi_pole/Phi_zero*Sig_zero ...
										) / Phi_zero;
									0,Phi_pole/Phi_zero];
							end
						end

					end
				end
			end

			[num,den] = tfdata(sys_R);
		    for jj=1:handles.n_i,
		        for ii=1:handles.n_o,
		        	n_c = size(TBa{ii+(jj-1)*n_o},2);

		            num{ii,jj}(end+1-[n_c:-1:1]) = ...
						(TBa{ii+(jj-1)*n_o}*...
							den{ii,jj}(end+1-[n_c:-1:1]).').';
				end
			end
			for jj=1:handles.n_i,
				for ii=1:handles.n_o,
                    sys_F = zpk(zF{ii+(jj-1)*n_o},pF, ...
    					1);
                    sys_est(ii,jj) = minreal( ...
                    	tf(num{ii,jj},den{ii,jj}) * sys_F);
                end
			end
		end

	%   Scaling
    for jj=1:n_i,
	    for ii=1:n_o,
			H_est = freqresp(sys_est(ii,jj), ...
				handles.freq(ind_f_eval)*2*pi);
			H_exp  = handles.H_exp(ind_f_eval,ii+(jj-1)*n_o);

%			sys_imsi(ii,jj) = sys_imsi(ii,jj)*abs(h_ref)/abs(h_imsi);
			sys_est(ii,jj) = sys_est(ii,jj)*H_exp/H_est;
	    end
	end

	handles.ind_sys = length(handles.sys_pool) + 1;
	handles.sys_pool{handles.ind_sys} = sys_est;
	str_ind_sys = [num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

	% Update handles structure
	guidata(hObject, handles);


% --- Executes on button press in pushbutton_Delete_Sys.
function pushbutton_Delete_Sys_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Delete_Sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.ind_sys > 1,
        handles.sys_pool(handles.ind_sys) = [];
    end
	handles.ind_sys = max([handles.ind_sys - 1,1]);

	str_ind_sys = [num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

	% Update handles structure
	guidata(hObject, handles);


% --- Executes on button press in pushbutton_Prev_Sys.
function pushbutton_Prev_Sys_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Prev_Sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.ind_sys = mod(handles.ind_sys - 2,length(handles.sys_pool)) + 1;

	str_ind_sys = [ ...
		num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in pushbutton_Next_Sys.
function pushbutton_Next_Sys_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Next_Sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.ind_sys = mod(handles.ind_sys,length(handles.sys_pool)) + 1;

	str_ind_sys = [ ...
		num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

    % Update handles structure
    guidata(hObject, handles);


function draw_n(handles)

    axes(handles.axes_mag) % Select the proper axes
	cla
	axes(handles.axes_pha) % Select the proper axes
	cla


function draw_frf(handles)

	str_ind_sys = [ ...
		num2str(handles.ind_f), ...
		'/', ...
		num2str(handles.n_io)];
	set(handles.text_ind_f, 'String', str_ind_sys);

	axes(handles.axes_pha) % Select the proper axes
   	plot(handles.freq,my_angle(handles.H_exp(:,handles.ind_f)), ...
		'Color',[0,0,0.75]);
	x_p = xlim;
    y_p = ylim;
    cla

	try,
		freq_est = [x_p(1):diff(x_p)/2^10:x_p(2)].';
		H_est_p = freqresp( ...
			handles.sys_pool{handles.ind_sys},freq_est*2*pi);
		H_est = [];
		for ii=1:handles.n_i,
		for jj=1:handles.n_o,
			H_est = [H_est, (squeeze(H_est_p(jj,ii,:)))];
		end
		end

	    axes(handles.axes_mag) % Select the proper axes
        set(gca,'XLimMode','auto');
        set(gca,'YLimMode','auto');
		hold on
%		semilogy(handles.freq,abs(H_mdl(:,handles.ind_f)), ...
%			'Color',[1.000 0.502 0.251]);
		semilogy(freq_est,abs(H_est(:,handles.ind_f)),'r');
		axes(handles.axes_pha) % Select the proper axes
		hold on
		grid off
%		plot(handles.freq,my_angle(H_mdl(:,handles.ind_f)), ...
%			'Color',[1.000 0.502 0.251]);
		plot(freq_est,my_angle(H_est(:,handles.ind_f)),'r');
	end

    axes(handles.axes_mag) % Select the proper axes
	semilogy(handles.freq,abs(handles.H_exp(:,handles.ind_f)), ...
		'Color',[0,0,0.75]);
	xlabel('Frequency (Hz)')
	ylabel('Magnitude')
	hold off
	grid on
	set(gca,'YMinorGrid','off')

	axes(handles.axes_pha) % Select the proper axes
	x_p = [min(handles.freq),max(handles.freq)];
    hold on
	grid off
	plot(x_p,[0,0],':k');
	plot(x_p,[-180,-180],':k');
	plot(x_p,[180,180],':k');
	plot(x_p,[-90,-90],':k');
	plot(x_p,[90,90],':k');
	plot(x_p,[-45,-45],':k');
	plot(x_p,[45,45],':k');
	plot(x_p,[-45,-45]-90,':k');
	plot(x_p,[45,45]+90,':k');
	set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
%	set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180])
	set(gca,'XGrid','on')
	plot(handles.freq,my_angle(handles.H_exp(:,handles.ind_f)), ...
		'Color',[0,0,0.75]);
	xlabel('Frequency (Hz)')
	ylabel('Phase (^o)')
    ylim(y_p);
	hold off


function pz_window(handles)

	sys_est = handles.sys_pool{handles.ind_sys};
	sys_R = sys_est;
	sys_F = sys_est;

	pF = handles.fixed_pole;
	zF = handles.fixed_zero;


%	if get(handles.togglebutton_Fixed_PZ,'value') > 0,
		for jj=1:handles.n_i,
			for ii=1:handles.n_o,
				sys_F(ii,jj) = zpk(zF{ii+(jj-1)*handles.n_o},pF, ...
					1);
				sys_R(ii,jj) = minreal(sys_est(ii,jj) / sys_F(ii,jj));
            end
		end
%	end

	[z,p,k] = zpkdata(sys_R);

	p = p{handles.ind_f};
	z = z{handles.ind_f};
	k = k(handles.ind_f);

	p = sort(p);
	z = sort(z);

		[zF,pF,kF] = zpkdata(sys_F);

		pF = pF{handles.ind_f};
		zF = zF{handles.ind_f};
		kF = kF(handles.ind_f);

		pF = sort(pF);
		zF = sort(zF);

	str1 = {};
	str2 = {};
	for ii=1:length(p),
%	    str1{ii} = sprintf('% 9.2g',real(p(ii)));
%	    str2{ii} = sprintf('% 9.2g',imag(p(ii)));
	    str1{ii} = num2str(real(p(ii)));
	    str2{ii} = num2str(imag(p(ii)));
	end

		for ii=1:length(pF),
	%	    str1{ii} = sprintf('% 9.2g',real(p(ii)));
	%	    str2{ii} = sprintf('% 9.2g',imag(p(ii)));
		    str1{ii+length(p)} = ['f',num2str(real(pF(ii)))];
		    str2{ii+length(p)} = ['f',num2str(imag(pF(ii)))];
		end

	set(handles.edit_Pole_R,'String',...
    	str1);
	set(handles.edit_Pole_I,'String',...
    	str2);

	str1 = {};
	str2 = {};
	for ii=1:length(z),
%	    str1{ii} = sprintf('% 9.2g',real(z(ii)));
%	    str2{ii} = sprintf('% 9.2g',imag(z(ii)));
	    str1{ii} = num2str(real(z(ii)));
	    str2{ii} = num2str(imag(z(ii)));
    end

		for ii=1:length(zF),
	%	    str1{ii} = sprintf('% 9.2g',real(p(ii)));
	%	    str2{ii} = sprintf('% 9.2g',imag(p(ii)));
		    str1{ii+length(z)} = ['f',num2str(real(zF(ii)))];
		    str2{ii+length(z)} = ['f',num2str(imag(zF(ii)))];
		end

	set(handles.edit_Zero_R,'String',...
    	str1);
	set(handles.edit_Zero_I,'String',...
    	str2);

	set(handles.edit_n_Poles,'String',...
    	num2str(length(p)+length(pF)));
	set(handles.edit_n_Zeros,'String',...
    	num2str(length(z)+length(zF)));

	draw_pz(handles);


function draw_pz(handles)

	[z,p,k] = zpkdata(handles.sys_pool{handles.ind_sys});

	p = p{handles.ind_f};
	z = z{handles.ind_f};
	k = k(handles.ind_f);

	p = sort(p);
	z = sort(z);

%	wD = imag(p);
%	xiN = real(p) ./ wN;

	fp = abs(p)/(2*pi);
	ind = isfinite(fp) & fp < max(handles.freq) & (imag(p) ~= 0);
	fp = fp(ind);

	fz = abs(z)/(2*pi);
	ind = isfinite(fz) & fz < max(handles.freq) & (imag(z) ~= 0);
	fz = fz(ind);

    axes(handles.axes_mag) % Select the proper axes
	y_p = get(handles.axes_mag,'YLim');
%	cla
	hold on
	for ii=1:length(fp),
	semilogy(fp(ii)*[1,1],y_p,':','Color',[0.5 0 0]);
	end
	for ii=1:length(fz),
	semilogy(fz(ii)*[1,1],y_p,':','Color',[0 0 0.5]);
	end
	hold off

	axes(handles.axes_pha) % Select the proper axes
	y_p = get(handles.axes_pha,'YLim');
%	cla
	hold on
	for ii=1:length(fp),
	plot(fp(ii)*[1,1],y_p,':','Color',[0.5 0 0]);
	end
	for ii=1:length(fz),
	plot(fz(ii)*[1,1],y_p,':','Color',[0 0 0.5]);
	end
	hold off


function edit_f_eval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f_eval as text
%        str2double(get(hObject,'String')) returns contents of edit_f_eval as a double

	handles.f_eval = str2double(get(hObject,'String'));

	%   Scaling
    [z,p,k] = zpkdata(handles.sys_pool{handles.ind_sys});
    ind_f_eval = sum(handles.freq <= handles.f_eval);
    set(handles.edit_f_eval, 'String', num2str(handles.f_eval));

    for ii=1:handles.n_io,
    sys_imsi = zpk(z{ii},p{ii},1);
    h_imsi = freqresp(sys_imsi,handles.freq(ind_f_eval)*2*pi);
    h_ref  = handles.H_exp(ind_f_eval,ii);

%	k(ii) = abs(h_ref)/abs(h_imsi);
    k(ii) = h_ref/h_imsi;
    end

	%handles.ind_sys = handles.ind_sys;
	handles.sys_pool{handles.ind_sys} = zpk(z,p,k);

	% Update handles structure
	guidata(hObject, handles);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end


% --- Executes during object creation, after setting all properties.
function edit_f_eval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_f_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton_Phase_Shift.
function pushbutton_Phase_Shift_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Phase_Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	[z,p,k] = zpkdata(handles.sys_pool{handles.ind_sys});
	k(handles.ind_f) = -k(handles.ind_f);

	%handles.ind_sys = handles.ind_sys;
	handles.sys_pool{handles.ind_sys} = zpk(z,p,k);

	% Update handles structure
	guidata(hObject, handles);

	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end


% --- Executes on button press in pushbutton_Data_Open.
function pushbutton_Data_Open_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Data_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	data_input = datain;

    handles.freq = data_input.freq;
    handles.H_exp = data_input.H_exp;
    handles.W_H = data_input.W_H;
    handles.W_H_inv = data_input.W_H_inv;

    handles.n_o = data_input.n_o;
    handles.n_i = data_input.n_i;
    handles.n_io = data_input.n_i*data_input.n_o;

    handles.n_pole = data_input.n_pole;
    handles.n_zero = data_input.n_zero*ones(handles.n_io,1);

    handles.ind_SR_o = data_input.ind_SR_o;
    handles.ind_SR_i = data_input.ind_SR_i;

	handles = data_pre(handles,0);

    guidata(hObject, handles);


% --- Executes on selection change in popupmenu_MDM_L.
function popupmenu_MDM_L_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_MDM_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_MDM_L contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_MDM_L

% load_listbox(pwd, handles);
% if strcmp(get(handles.fig_mfdid,'SelectionType'),'open') % If double click
    index_selected = get(handles.popupmenu_MDM_L,'Value');
    file_list = get(handles.popupmenu_MDM_L,'String');
    filename = file_list{index_selected}; % Item selected in list box
    try,
        cd (filename)
        load_listbox(pwd,handles) % Load list box with new directory
    catch,
        [path,name,ext] = fileparts(filename);
        switch ext
        case '.mat'
            load(filename) % Open FIG-file with guide

            handles.ind_sys = handles.ind_sys + 1;
            handles.sys_pool{handles.ind_sys} = sys_mfdid;

			str_ind_sys = [ ...
				num2str(handles.ind_sys), ...
				'/', ...
				num2str(length(handles.sys_pool))];
			set(handles.text_ind_sys, 'String', str_ind_sys);

			n_error_norm(handles);

            draw_n(handles);
			draw_frf(handles);
    		try,
			pz_window(handles);
			end

			% Update handles structure
			guidata(hObject, handles);

		otherwise
        	try
        		open(filename) % Use open for other file types
            catch
            	errordlg(lasterr,'File Type Error','modal')
			end
        end
    end
% end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_MDM_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_MDM_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
	dir_struct = dir(pwd);
	[sorted_names,sorted_index] = sortrows({dir_struct.name}');
	set(hObject,'String',sorted_names,...
	    'Value',1)


function edit_MDM_S_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MDM_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MDM_S as text
%        str2double(get(hObject,'String')) returns contents of edit_MDM_S as a double

assignin('base','sys_mfdid',handles.sys_pool{handles.ind_sys});
[zz,pp,kk] = zpkdata(handles.sys_pool{handles.ind_sys});
assignin('base','z',zz);
assignin('base','p',pp);
assignin('base','k',kk);
evalin('base',['save ', get(hObject,'String'), ' sys_mfdid ','z ','p ','k']);


% --- Executes during object creation, after setting all properties.
function edit_MDM_S_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MDM_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton_Fr_Add.
function pushbutton_Fr_Add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Fr_Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ind_freq = find(handles.ind_freq);
d_ind_freq = diff(ind_freq);
ind0 = [ind_freq(1);ind_freq(find(d_ind_freq-1)+1)];
ind1 = [ind_freq(find(d_ind_freq-1));ind_freq(end)];

y01 = get(handles.axes_pha,'YLim');
axes(handles.axes_pha) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        [-200,200],'k:')
end
ind_tick = [-180,-135,-90,-45,0,45,90,135,180];
for ii=1:length(ind_tick),
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

y01 = get(handles.axes_mag,'YLim');
axes(handles.axes_mag) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        get(gca,'Ylim'),'k:')
end
ind_tick = get(gca,'YTick');
for ii=2:length(ind_tick)-1,
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

%
	imsi = ginput(2);

	f0 = imsi(1,1);
	ff = imsi(2,1);	clear imsi

	handles.ind_freq((handles.freq >= f0) & (handles.freq <= ff)) ...
        = 1;

%
	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

ind_freq = find(handles.ind_freq);
d_ind_freq = diff(ind_freq);
ind0 = [ind_freq(1);ind_freq(find(d_ind_freq-1)+1)];
ind1 = [ind_freq(find(d_ind_freq-1));ind_freq(end)];

y01 = get(handles.axes_pha,'YLim');
axes(handles.axes_pha) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        [-200,200],'k:')
end
ind_tick = [-180,-135,-90,-45,0,45,90,135,180];
for ii=1:length(ind_tick),
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

y01 = get(handles.axes_mag,'YLim');
axes(handles.axes_mag) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        get(gca,'Ylim'),'k:')
end
ind_tick = get(gca,'YTick');
for ii=2:length(ind_tick)-1,
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

guidata(hObject, handles);


% --- Executes on button press in pushbutton_Fr_Del.
function pushbutton_Fr_Del_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Fr_Del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ind_freq = find(handles.ind_freq);
d_ind_freq = diff(ind_freq);
ind0 = [ind_freq(1);ind_freq(find(d_ind_freq-1)+1)];
ind1 = [ind_freq(find(d_ind_freq-1));ind_freq(end)];

y01 = get(handles.axes_pha,'YLim');
axes(handles.axes_pha) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        [-200,200],'k:')
end
ind_tick = [-180,-135,-90,-45,0,45,90,135,180];
for ii=1:length(ind_tick),
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

y01 = get(handles.axes_mag,'YLim');
axes(handles.axes_mag) % Select the proper axes
hold on
grid off
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
		[y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        get(gca,'Ylim'),'k:')
end
ind_tick = get(gca,'YTick');
for ii=2:length(ind_tick)-1,
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

%
	imsi = ginput(2);

	f0 = imsi(1,1);
	ff = imsi(2,1);	clear imsi

	handles.ind_freq((handles.freq >= f0) & (handles.freq <= ff)) ...
        = 0;

%
	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

ind_freq = find(handles.ind_freq);
d_ind_freq = diff(ind_freq);
ind0 = [ind_freq(1);ind_freq(find(d_ind_freq-1)+1)];
ind1 = [ind_freq(find(d_ind_freq-1));ind_freq(end)];

y01 = get(handles.axes_pha,'YLim');
axes(handles.axes_pha) % Select the proper axes
hold on
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
        [y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
hold off

y01 = get(handles.axes_mag,'YLim');
axes(handles.axes_mag) % Select the proper axes
hold on
for ii=1:length(ind0),
%	h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
%		[y01(1),y01(2),y01(2),y01(1)],[1 1 0.75]);
    h_fill = fill(handles.freq([ind0(ii),ind0(ii),ind1(ii),ind1(ii)]), ...
		[y01(1),y01(2),y01(2),y01(1)],[ 0.875 1.000 0.875 ]);
    set(h_fill,'FaceAlpha',0.5)
end
ind_tick = get(gca,'XTick');
for ii=2:length(ind_tick)-1,
    plot([ind_tick(ii),ind_tick(ii)], ...
        get(gca,'Ylim'),'k:')
end
ind_tick = get(gca,'YTick');
for ii=2:length(ind_tick)-1,
    plot(get(gca,'Xlim'), ...
        [ind_tick(ii),ind_tick(ii)],'k:')
end
hold off

guidata(hObject, handles);


% --- Executes on button press in pushbutton_Modal.
function pushbutton_Modal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Modal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	draw_n(handles);
	draw_frf(handles);
    try,
	pz_window(handles);
	end

	[z,p,k] = zpkdata(handles.sys_pool{handles.ind_sys});

	p = p{handles.ind_f};
	z = z{handles.ind_f};
	k = k(handles.ind_f);

    p = sort(p);
    z = sort(z);

    p = p(imag(p) ~= 0);
    z = z(imag(z) ~= 0);

    axes(handles.axes_mag) % Select the proper axes
%   xy = get(hObject,'CurrentPoint');
%   xy = get(GCB0,'CurrentPoint');
    [xp,yp] = ginput(1);

    [mp,mpi] = min(abs(abs(p)/(2*pi) - xp*ones(size(p))));
    [mz,mzi] = min(abs(abs(z)/(2*pi) - xp*ones(size(z))));

    if  ~isempty(mp) & ...
        (all(mp < mz) | isempty(mz)),
        ind_pz = 'pole';
        fN = abs(p(mpi))/(2*pi);
        xiN = -real(p(mpi)) ./ (fN*2*pi);
    elseif ~isempty(mz),
        ind_pz = 'zero';
        fN = abs(z(mzi))/(2*pi);
        xiN = -real(z(mzi)) ./ (fN*2*pi);
    end

    try,
    str = sprintf('%s\nfN : %g\nxiN : %g',ind_pz,fN,xiN);
    xl = get(handles.axes_mag,'XLim');
    xp = min([xl(2),fN]);
    yp = freqresp(handles.sys_pool{handles.ind_sys},fN*2*pi);
    yp = abs(yp(handles.ind_f));
    hold on
    plot(xp,yp,'o')
    ht = text(xp+xl(2)*0.02,yp,str);
    set(ht,'Color','m','FontAngle','oblique','FontWeight','Bold')
%         'BackgroundColor','w');
    alpha(0.5)
    end

if 0,
    figure
    uicontrol('Style', 'Text', ...
         'String', 'Pole : Frequency', ...
         'Position', [20,390,100,50])
    uicontrol('Style', 'Text', ...
         'String', 'Pole : Damping Ratio', ...
         'Position', [120,390,100,50])
    uicontrol('Style', 'Text', ...
         'String', 'Zero : Frequency', ...
         'Position', [220,390,100,50])
    uicontrol('Style', 'Text', ...
         'String', 'Zero : Damping Ratio', ...
         'Position', [320,390,100,50])

    str_pw = sprintf('\n%g',abs(p)/(2*pi));
    str_px = sprintf('\n%g',-real(p)./abs(p));
    str_zw = sprintf('\n%g',abs(z)/(2*pi));
    str_zx = sprintf('\n%g',-real(z)./abs(z));
end

fprintf(1,'\nPole\nwN\txi');
for ii=1:size(p),
fprintf(1,'\n%g\t%g',abs(p(ii))/(2*pi),-real(p(ii))./abs(p(ii)));
end
fprintf(1,'\n\nZero\nwN\txi');
for ii=1:size(z),
fprintf(1,'\n%g\t%g',abs(z(ii))/(2*pi),-real(z(ii))./abs(z(ii)));
end

if 0,
    uicontrol('Style', 'ListBox', ...
         'String', str_pw, ...
         'Position', [20,20,100,350])
    uicontrol('Style', 'ListBox', ...
         'String', str_px, ...
         'Position', [130,20,100,350])
    uicontrol('Style', 'ListBox', ...
         'String', str_zw, ...
         'Position', [240,20,100,350])
    uicontrol('Style', 'ListBox', ...
         'String', str_zx, ...
         'Position', [350,20,100,350])
end


% --- Executes on button press in togglebutton_Inv.
function togglebutton_Inv_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_Inv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_Inv

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
    set(hObject,'FontWeight','bold');
    set(hObject,'ForegroundColor',[1,1,0]);
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
    set(hObject,'FontWeight','normal');
    set(hObject,'ForegroundColor',[0,0,0]);
end


% --- Executes on button press in pushbutton_PZ_plot.
function pushbutton_PZ_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PZ_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	[zN,pN,kN] = zpkdata(handles.sys_pool{handles.ind_sys});

h_f_pz = figure;
plot(real(pN{handles.ind_f}),imag(pN{handles.ind_f}), ...
    'x','MarkerSize',6,'Color',[0,0,0]);
hold on
plot(real(zN{handles.ind_f}),imag(zN{handles.ind_f}), ...
    'o','MarkerSize',6,'Color',[1,0,0]);
hold off


% --- Executes on button press in pushbutton_GN.
function pushbutton_GN_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%
%	PART I.		Pre Process
%

	%
	%	Verification of Setting
	%
	ind_verb = 0;

	handles.ind_ck = 0;
	if handles.ind_ck == 0,
		warning off
		echo off
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	%
	%	Basic Parameter Setting
	%
	im = sqrt(-1);

	freq     = handles.freq(handles.ind_freq);
	H_exp    = handles.H_exp(handles.ind_freq,:);
	H_exp_inv= 1./H_exp;
	n_o      = handles.n_o;
	n_i      = handles.n_i;
	n_io     = handles.n_io;
	ind_SR_o = handles.ind_SR_o;
	ind_SR_i = handles.ind_SR_i;
	W_H      = handles.W_H(handles.ind_freq,:);
%	W_H      = ones(size(handles.W_H(handles.ind_freq,:)));
	W_H_inv  = handles.W_H_inv(handles.ind_freq,:);
%	W_H_inv  = ones(size(handles.W_H_inv(handles.ind_freq,:)));
	n_pole	 = handles.n_pole;
	n_zero	 = handles.n_zero;
	n_pole_F = 0;
	n_zero_F = zeros(n_io,1);
    n_pole_R = n_pole - n_pole_F;
    n_zero_R = n_zero - n_zero_F;
	pF		 = handles.fixed_pole;
	zF		 = handles.fixed_zero;

	s = im*2*pi*freq;
    n_freq = length(freq);

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(H_exp))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(W_H))))])
		pause
	end

	%
	%	Option : Real Flag
	%
	ind_Real = 0;
    if get(handles.togglebutton_Real,'Value') > 0,
	    ind_Real = 1;
	end

	%
	%	Option : Fixed
	%
	ind_Fixed = 0;
    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
	    ind_Fixed = 1;

		for ii=1:length(pF),
			H_exp = H_exp ...
				.* ...
				((s - pF(ii)) ...
					*ones(1,n_io));
			W_H = W_H ...
				./ ...
				((s - pF(ii)) ...
					*ones(1,n_io));
		end
		for jj=1:n_io,
		for ii=1:length(zF{jj}),
			H_exp(:,jj) = H_exp(:,jj) ...
				./ (s - zF{jj}(ii));
			W_H(:,jj) = W_H(:,jj) ...
				.* (s - zF{jj}(ii));
		end
		end

	    ind_sc = ~isfinite(H_exp);
	    H_exp(ind_sc) = 0;
		W_H(ind_sc) = 0;

	    n_pole_F = length(pF);
		for ii=1:n_io,
			n_zero_F(ii) = length(zF{ii});
	    end
	    n_pole_R = n_pole - n_pole_F;
	    n_zero_R = n_zero - n_zero_F;
    end

	%
	%	Option : Structural Relationships
	%
	ind_SR = 0;
    TBa = cell(n_io,1);
    if get(handles.togglebutton_SR,'Value') > 0,
	    ind_SR = 1;

	    for ii=1:n_i,
		    for jj=1:n_o,
		    	switch ind_SR_i(ii),
		    		case 0,
		    			TBa{jj+(ii-1)*n_o} = [];
		    		case 1,
		    			switch ind_SR_o(jj),
		    				case {0,1,4},
		    					TBa{jj+(ii-1)*n_o} = [];
		    				case {2,5},
		    					TBa{jj+(ii-1)*n_o} = [0];
		    				case {3},
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
                            case {6},
                                TBa{jj+(ii-1)*n_o} = zeros(4,4);
		    			end
		    		case 2,
		    			switch ind_SR_o(jj),
		    				case {0,1},
		    					TBa{jj+(ii-1)*n_o} = [];
		    				case {2},
		    					TBa{jj+(ii-1)*n_o} = [0];
		    				case {3},
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
		    				case {4,5,6},
		    					TBa{jj+(ii-1)*n_o} = eye(2,2);
		    			end
		    	end
			end
		end

		if ind_Fixed == 1,
			Phi_pole = prod(-pF);
			Sig_pole = 0;
			for ii=1:length(pF),
				fixed_pole = pF;
				fixed_pole(ii) = [];
				Sig_pole = Sig_pole + prod(-fixed_pole);
			end

		    for ii=1:n_i,
			    for jj=1:n_o,

					if ((ind_SR_i(ii) == 2) & (ind_SR_o(jj) > 3)),
				    	fixed_zero = zF{jj+(ii-1)*n_o};
						Phi_zero = ...
							prod(-fixed_zero);
						Sig_zero = 0;

						if Phi_zero == 0,
	    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
						else,
							for kk = 1:length(fixed_zero),
								fixed_zero_kk  = fixed_zero;
								fixed_zero_kk(kk) = [];
								Sig_zero = Sig_zero ...
									+ prod(-fixed_zero_kk);
							end
	    					TBa{jj+(ii-1)*n_o} = ...
								[Phi_pole/Phi_zero, ...
									( ...
										Sig_pole ...
										-Phi_pole/Phi_zero*Sig_zero ...
									) / Phi_zero;
								0,Phi_pole/Phi_zero];
						end
					end

				end
			end
		end
    end

	%
	%	Option : Inverse Model
	%
    ind_Inv = 0;
    if get(handles.togglebutton_Inv,'Value') > 0,
	    ind_Inv = 1;

		if ind_Fixed == 1,
			for ii=1:length(pF),
				H_exp_inv = H_exp_inv ...
					./ ...
					( ...
						(s - pF(ii)) ...
						*ones(1,n_io) ...
					);

				W_H_inv = W_H_inv ...
					.* ...
					( ...
						(s - pF(ii)) ...
						*ones(1,n_io) ...
					);
			end
			for jj=1:n_io,
				for ii=1:length(zF{jj}),
					W_H_inv(:,jj) = W_H_inv(:,jj) ...
						./ (s - zF{jj}(ii));
					H_exp_inv(:,jj) = H_exp_inv(:,jj) ...
						.* (s - zF{jj}(ii));
				end
			end
		end

		ind_sc = ~isfinite(H_exp_inv);
		H_exp_inv(ind_sc) = 0;
		W_H_inv(ind_sc) = 0;
    end


%
%	PART II.	Main Process
%

	%
	%	Stacked Matrices for the System
	%
	H_exp = reshape(H_exp,n_freq*n_io,1);
	W_H	  = reshape(W_H  ,n_freq*n_io,1);

	if ind_Inv == 1,
		H_exp_inv = reshape(H_exp_inv,n_freq*n_io,1);
		W_H_inv	  = reshape(W_H_inv  ,n_freq*n_io,1);
	end

	sE = [];
	for ii=1:n_io,
		sE = [sE; ...
			s];
	end

	n_max = max(n_pole_R, max(n_zero_R) + 1);
	Phi = ones(n_freq,n_max);
	for ii = 2:n_max,
		Phi(:,ii) = s.^(ii-1);
	end

	%
	%	ARMA Parameters
	%
		%
		%	Phi_a
		%
		Phi_a = zeros(n_freq*n_io,n_pole_R);
		for ii = 1:n_io,
			Phi_a([1:n_freq]+(ii-1)*n_freq,:) = ...
				Phi(:,[n_pole_R:-1:1]);
		end

		%
		%	Phi_B
		%
		Phi_B = zeros(n_freq*n_io,sum(n_zero_R)+n_io);
		for ii=1:n_io,
			ind_B = (n_zero_R(ii)+1):-1:1;
			Phi_B([1:n_freq]+(ii-1)*n_freq, ...
				[1:(n_zero_R(ii)+1)] + ...
				ii-1+sum(n_zero_R(1:ii-1))) = ...
				Phi(:,ind_B);
		end

		if ind_SR == 1,
		%
		%	Phi_a_c & Phi_B_c
		%
		Phi_a_c = zeros(n_freq*n_io,n_pole_R);
		Phi_B_c = Phi_B;
		for ii=n_io:-1:1,
			n_c = size(TBa{ii},2);
			Phi_a_c([1:n_freq]+(ii-1)*n_freq, ...
				(n_pole_R+1-[n_c:-1:1])) = ...
				Phi(:,[n_c:-1:1]) * TBa{ii};
			Phi_B_c(:,n_zero_R(ii)+2+ ...
				sum(n_zero_R(1:ii-1))+(ii-1)- ...
				[n_c:-1:1]) = [];
		end
		end

	%
	%	Initial Parameter
	%
	sys_R = handles.sys_pool{handles.ind_sys};
	sys_F = sys_R;

	if ind_Fixed == 1,
		for jj=1:handles.n_i,
			for ii=1:handles.n_o,
				sys_F(ii,jj) = zpk(zF{ii+(jj-1)*n_o},pF, ...
					1);
				sys_R(ii,jj) = minreal(sys_R(ii,jj) ...
					/ sys_F(ii,jj));
            end
		end
	end
	[num,den] = tfdata(sys_R);

	theta_a = den{1,1}.';

	ind_srt = length(theta_a) - n_pole_R;
	if ind_srt > 1,
		if handles.ind_ck == 1,
	    disp('check coef of pole')
	    disp(theta_a(1:ind_srt-1));
	    end
	    theta_a = theta_a(ind_srt:end);
	end

    theta_B = [];
	for jj=1:n_i,
		for ii=1:n_o,
            theta_B_imsi = num{ii,jj}.' / theta_a(1);

            ind_srt = length(theta_B_imsi) ...
                - n_zero_R(ii+(jj-1)*n_o);
            if ind_srt > 1,
            	if handles.ind_ck == 1,
                disp('check coef of zero')
                disp(theta_B_imsi(1:ind_srt-1));
                end
                theta_B_imsi = theta_B_imsi(ind_srt:end);
            end

            if ind_SR == 0,
	            theta_B = [theta_B; ...
					theta_B_imsi];
			else,
				n_c = size(TBa{ii+(jj-1)*n_o},2);
	            theta_B = [theta_B; ...
					theta_B_imsi([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c])];
			end
		end
	end
	theta_a = theta_a / theta_a(1);
	theta_a(1) = [];

	theta = [theta_a;theta_B];

	%
	%	Initial Estimate
	%
	if ind_SR == 0,
		if n_pole_R == 0,
			H_est = (Phi_B * theta_B);
		else,
			H_est = (Phi_B * theta_B) ...
				./ (sE.^n_pole_R + Phi_a * theta_a);
		end
	else,
		if n_pole_R == 0,
			H_est = (Phi_B_c * theta_B ...
				+ Phi_a_c * theta_a);
		else,
			H_est = (Phi_B_c * theta_B ...
				+ Phi_a_c * theta_a) ...
				./ (sE.^n_pole_R + Phi_a * theta_a);
		end
	end

	eV = (H_exp-H_est).*W_H;
	if ind_Inv == 1,
		H_est_inv = 1 ./ H_est;
		eV = [eV;
			(H_exp_inv-H_est_inv).*W_H_inv];
	end
	ind_ck = ~isfinite(eV);
%	ind_ck = logical(sign(sum(ind_ck.')));
	eV(ind_ck) = 0;

	eN = eV'*eV;
	eN_pool = [eN];

	eN_min = eN;

	if (ind_verb),
%	    clc
	    disp(['INITIAL ESTIMATE'])
		disp(['Current fit: ' num2str(eN)])
%		disp(['par-vector'])
%		disp(theta)
	end

	h_f_gn = figure;
	semilogy(0,eN_pool(1),'or')
	xlabel('Iteration #')
	ylabel('Error Norm')
	title('Gauss-Newton Method')
	hold on

	%
	%	Convergence Conditions
	%
	tol = 0.01;

	tol_r = tol *ones(size(theta));
	tol_a = tol *ones(size(theta));
	tol_b = zeros(size(theta));
	tol_n = 30;

	%
	%	Iteration
	%
	st = 0;
	for ind_iter=1:tol_n,

		%	Update
    	theta_o = theta;
    	eN_o = eN;

		if n_pole_R > 0,
			den_R = (sE.^n_pole_R + Phi_a * theta_a);
		else,
			den_R = 1;
		end

	    %	Gradient
		if ind_SR == 0,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B,2);

			Del1_eV = ((H_est ./ den_R) ...
				*ones(1,ind_K1)) .* Phi_a;

			Del2_eV = - Phi_B ./ ...
				(den_R*ones(1,ind_K2));

			Del_eV = (W_H*ones(1,ind_K1 + ind_K2)) ...
				.*[Del1_eV Del2_eV];
		else,
			ind_K1 = size(Phi_a,2);
			ind_K2 = size(Phi_B_c,2);

			Del1_eV = ((H_est ./ den_R) ...
				*ones(1,ind_K1)) .* Phi_a ...
				- Phi_a_c ./ ...
                (den_R*ones(1,ind_K1));

			Del2_eV = - Phi_B_c ./ ...
				(den_R*ones(1,ind_K2));

			Del_eV = (W_H*ones(1,ind_K1 + ind_K2)) ...
				.*[Del1_eV Del2_eV];
		end

		if ind_Inv == 1,
			if ind_SR == 0,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B,2);

				Del1_eV_inv = -Phi_a ...
					./ (Phi_B * theta_B *ones(1,ind_K1));

				Del2_eV_inv = ...
					(( ...
						H_est_inv ...
						./ ...
						(Phi_B * theta_B) ...
					)*ones(1,ind_K2)) ...
					.* Phi_B;

				Del_eV_inv = (W_H_inv*ones(1,ind_K1 + ind_K2)) ...
					.*[Del1_eV_inv Del2_eV_inv];
			else,
				ind_K1 = size(Phi_a,2);
				ind_K2 = size(Phi_B_c,2);

				Del1_eV_inv = -Phi_a ...
					./ ((Phi_B_c * theta_B + Phi_a_c * theta_a) ...
						*ones(1,ind_K1)) ...
					+ ( ...
						(H_est_inv ...
						./ ...
						(Phi_B_c * theta_B + Phi_a_c * theta_a)) ...
						*ones(1,ind_K1)) ...
						.* Phi_a_c;

				Del2_eV_inv = ( ...
					(H_est_inv ...
					./ ...
					(Phi_B_c * theta_B +  Phi_a_c * theta_a)) ...
					* ones(1,ind_K2)) ...
					.* Phi_B_c;

				Del_eV_inv = (W_H_inv*ones(1,ind_K1 + ind_K2)) ...
					.*[Del1_eV_inv Del2_eV_inv];
			end
			Del_eV = [Del_eV;
				Del_eV_inv];
		end

		ind_ck = ~isfinite(Del_eV);
		ind_ck = logical(sign(sum(ind_ck.')));
		Del_eV(ind_ck,:) = 0;

		%	Gauss-Newton Search Direction
		K = Del_eV'*Del_eV;
		f = Del_eV'*eV;
		if ind_Real
			K = real(K);
			f = real(f);
		end
		del_th = K\f;

	    %	Line Search
		alpha_GN = 1;
		for ind_iter2=1:tol_n,
			theta = theta_o - alpha_GN*del_th;
			if ind_iter2 == tol_n-1, theta = theta_o; end

			%	Parameters
			theta_a = [1; ...
				theta(1:n_pole_R)];
	    	theta_a = apolystab(theta_a,ind_Real);
		    theta_B = theta(n_pole_R+1:end);

			theta_a(1) = [];
			theta = [theta_a;theta_B];

			%	Error Norm Check
            if ind_SR == 0,
				H_est = (Phi_B * theta_B) ...
					./ (sE.^n_pole_R + Phi_a * theta_a);
			else,
				H_est = (Phi_B_c * theta_B ...
					+ Phi_a_c * theta_a) ...
					./ (sE.^n_pole_R + Phi_a * theta_a);
			end

			eV = (H_exp-H_est).*W_H;

			if ind_Inv == 1,
				H_est_inv = 1 ./ H_est;
				eV = [eV;
					(H_exp_inv-H_est_inv).*W_H_inv];
			end

			ind_ck = ~isfinite(eV);
%			ind_ck = logical(sign(sum(ind_ck.')));
			eV(ind_ck) = 0;

			eN = eV'*eV;
			eN_pool = [eN_pool;eN];

			if (ind_verb),
	            home, disp(int2str(ind_iter2))
	        end;

			figure(h_f_gn);
			semilogy(length(eN_pool)-1,eN_pool(end),'o')

            if eN < eN_min,
                theta_a_min = theta_a;
                theta_B_min = theta_B;
                eN_min = eN;
            end

	    	if eN < eN_o,
	    		break
	    	end

	        if ind_iter2==tol_n,st=1;end

			%	Shooting Again !!!
			alpha_GN = alpha_GN/2;

	        if ind_iter2 == ceil(tol_n/2),
	        	del_th = f/norm(K)*length(K);
	        	alpha_GN = 1;
	        end
	    end

		%	Convergence Check
	    if (ind_verb),
	        home
	        disp(['ITERATION # ' int2str(ind_iter)])
	        disp(['Current fit:  ', num2str(eN), ...
	        	'  Previous fit:  ' num2str(eN_o)])
%			disp(['Current par prev par   GN-dir'])
%			disp([theta theta_o del_th])
	        disp(['Norm of GN-vector: ' num2str(norm(del_th))])
	        if st==1,
	            disp(['No improvement of the criterion possible', ...
	            	' along the search direction.'])
	            disp('Iterations therefore terminated.')
	        end
	    end

		figure(h_f_gn);
		semilogy(length(eN_pool)-1,eN_pool(end),'ok')

		tol_b = max(tol_b,theta);
		if (abs((theta - theta_o)./theta_o) < tol_r) ...
			& (abs(theta - theta_o) < tol_a.*tol_b),
			break
		end

		if (st == 1) | (norm(del_th) < tol_a),
			break
		end
	end

	figure(h_f_gn);
	semilogy(length(eN_pool),eN_min(end),'or')
	hold off


%
%	PART III.	Post Process
%
	%
	%	Parameter Vectors
	%
	theta_a_R = [1;theta(1:n_pole_R)];
   	theta_a_R = apolystab(theta_a_R,ind_Real);
    theta_B_R = cell(n_io,1);
	if ind_SR == 0,
		for jj=1:n_i,
			for ii=1:n_o,
				theta_B_R{ii+(jj-1)*n_o} = ...
					theta([1:(n_zero_R(ii+(jj-1)*n_o)+1)] + ...
						sum(n_zero_R(1:(ii-1+(jj-1)*n_o))) + ...
						(ii-1+(jj-1)*n_o) + ...
						n_pole_R);
			end
		end
	else,
		ind_th = n_pole_R;
		for jj=1:n_i,
			for ii=1:n_o,
				n_c = size(TBa{ii+(jj-1)*n_o},2);
				theta_B_R{ii+(jj-1)*n_o} = ...
					[theta([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c] + ...
						ind_th);
						TBa{ii+(jj-1)*n_o}*...
							theta_a_R(n_pole_R+2-[n_c:-1:1]);
					];
				ind_th = ind_th + ...
					n_zero_R(ii+(jj-1)*n_o)+1-n_c;
			end
		end
	end

	%
	%	Rational Polynomial Model
	%
	for jj=1:n_i,
		for ii=1:n_o,
			num{ii,jj} = theta_B_R{ii+(jj-1)*n_o}.';
			den{ii,jj} = theta_a_R.';
		end
	end
	sys_R = tf(num,den);

	if ind_Fixed == 1,
		[z,p,k] = zpkdata(sys_R);
		for jj=1:n_i,
			for ii=1:n_o,
        		p{ii,jj} = [p{ii,jj}; ...
                    handles.fixed_pole];
				z{ii,jj} = [z{ii,jj}; ...
                        handles.fixed_zero{ii+(jj-1)*n_o}];
			end
		end
		sys_R = zpk(z,p,k);
	else,
		handles.fixed_pole = [];
		for ii = 1:handles.n_io,
			handles.fixed_zero{ii} = [];
		end
	end

	%
	%	Data Base & Graphic
	%
    handles.ind_sys = length(handles.sys_pool) + 1;
    handles.sys_pool{handles.ind_sys} = ...
 		sys_R;

	str_ind_sys = [ ...
		num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);
	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	    pz_window(handles);
	end

		%   Update handles structure
		guidata(hObject, handles);


% --- Executes on button press in pushbutton_LM.
function pushbutton_LM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	global eN_pool

%
%	PART I.		Pre Process
%

	%
	%	Verification of Setting
	%
	ind_verb = 0;

	handles.ind_ck = 0;
	if handles.ind_ck == 0,
		warning off
		echo off
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(handles.H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(handles.W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(handles.W_H))))])
		pause
	end

	%
	%	Basic Parameter Setting
	%
	im = sqrt(-1);

	freq     = handles.freq(handles.ind_freq);
	H_exp    = handles.H_exp(handles.ind_freq,:);
	H_exp_inv= 1./H_exp;
	n_o      = handles.n_o;
	n_i      = handles.n_i;
	n_io     = handles.n_io;
	ind_SR_o = handles.ind_SR_o;
	ind_SR_i = handles.ind_SR_i;
	W_H      = handles.W_H(handles.ind_freq,:);
%	W_H      = ones(size(handles.W_H(handles.ind_freq,:)));
	W_H_inv  = handles.W_H_inv(handles.ind_freq,:);
%	W_H_inv  = ones(size(handles.W_H_inv(handles.ind_freq,:)));
	n_pole	 = handles.n_pole;
	n_zero	 = handles.n_zero;
	n_pole_F = 0;
	n_zero_F = zeros(n_io,1);
    n_pole_R = n_pole - n_pole_F;
    n_zero_R = n_zero - n_zero_F;
	pF		 = handles.fixed_pole;
	zF		 = handles.fixed_zero;

	s = im*2*pi*freq;
    n_freq = length(freq);

	if handles.ind_ck > 0,
		figure(handles.ind_ck)
		waterfall(log(H_exp).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['H_{exp} ', ...
			num2str(sum(sum(~isfinite(H_exp))))])
		pause
	end

	if handles.ind_ck > 0,
		figure(handles.ind_ck+1)
		waterfall(log(W_H).')
		xlabel('Frequency Line')
		ylabel('Channel Line')
		title(['W_H ', ...
			num2str(sum(sum(~isfinite(W_H))))])
		pause
	end

	%
	%	Option : Real Flag
	%
	ind_Real = 0;
    if get(handles.togglebutton_Real,'Value') > 0,
	    ind_Real = 1;
	end

	%
	%	Option : Fixed
	%
	ind_Fixed = 0;
    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
	    ind_Fixed = 1;

		for ii=1:length(pF),
			H_exp = H_exp ...
				.* ...
				((s - pF(ii)) ...
					*ones(1,n_io));
			W_H = W_H ...
				./ ...
				((s - pF(ii)) ...
					*ones(1,n_io));
		end
		for jj=1:n_io,
		for ii=1:length(zF{jj}),
			H_exp(:,jj) = H_exp(:,jj) ...
				./ (s - zF{jj}(ii));
			W_H(:,jj) = W_H(:,jj) ...
				.* (s - zF{jj}(ii));
		end
		end

	    ind_sc = ~isfinite(H_exp);
	    H_exp(ind_sc) = 0;
		W_H(ind_sc) = 0;

	    n_pole_F = length(pF);
		for ii=1:n_io,
			n_zero_F(ii) = length(zF{ii});
	    end
	    n_pole_R = n_pole - n_pole_F;
	    n_zero_R = n_zero - n_zero_F;
    end

	%
	%	Option : Structural Relationships
	%
	ind_SR = 0;
    TBa = cell(n_io,1);
    if get(handles.togglebutton_SR,'Value') > 0,
	    ind_SR = 1;

	    for ii=1:n_i,
		    for jj=1:n_o,
		    	switch ind_SR_i(ii),
		    		case 0,
		    			TBa{jj+(ii-1)*n_o} = [];
		    		case 1,
		    			switch ind_SR_o(jj),
		    				case {0,1,4},
		    					TBa{jj+(ii-1)*n_o} = [];
		    				case {2,5},
		    					TBa{jj+(ii-1)*n_o} = [0];
		    				case {3},
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
                            case {6},
                                TBa{jj+(ii-1)*n_o} = zeros(4,4);
		    			end
		    		case 2,
		    			switch ind_SR_o(jj),
		    				case {0,1},
		    					TBa{jj+(ii-1)*n_o} = [];
		    				case {2},
		    					TBa{jj+(ii-1)*n_o} = [0];
		    				case {3},
		    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
		    				case {4,5,6},
		    					TBa{jj+(ii-1)*n_o} = eye(2,2);
		    			end
		    	end
			end
		end

		if ind_Fixed == 1,
			Phi_pole = prod(-pF);
			Sig_pole = 0;
			for ii=1:length(pF),
				fixed_pole = pF;
				fixed_pole(ii) = [];
				Sig_pole = Sig_pole + prod(-fixed_pole);
			end

		    for ii=1:n_i,
			    for jj=1:n_o,

					if ((ind_SR_i(ii) == 2) & (ind_SR_o(jj) > 3)),
				    	fixed_zero = zF{jj+(ii-1)*n_o};
						Phi_zero = ...
							prod(-fixed_zero);
						Sig_zero = 0;

						if Phi_zero == 0,
	    					TBa{jj+(ii-1)*n_o} = zeros(2,2);
						else,
							for kk = 1:length(fixed_zero),
								fixed_zero_kk  = fixed_zero;
								fixed_zero_kk(kk) = [];
								Sig_zero = Sig_zero ...
									+ prod(-fixed_zero_kk);
							end
	    					TBa{jj+(ii-1)*n_o} = ...
								[Phi_pole/Phi_zero, ...
									( ...
										Sig_pole ...
										-Phi_pole/Phi_zero*Sig_zero ...
									) / Phi_zero;
								0,Phi_pole/Phi_zero];
						end
					end

				end
			end
		end
    end

	%
	%	Option : Inverse Model
	%
    ind_Inv = 0;
    if get(handles.togglebutton_Inv,'Value') > 0,
	    ind_Inv = 1;

		if ind_Fixed == 1,
			for ii=1:length(pF),
				H_exp_inv = H_exp_inv ...
					./ ...
					( ...
						(s - pF(ii)) ...
						*ones(1,n_io) ...
					);

				W_H_inv = W_H_inv ...
					.* ...
					( ...
						(s - pF(ii)) ...
						*ones(1,n_io) ...
					);
			end
			for jj=1:n_io,
				for ii=1:length(zF{jj}),
					W_H_inv(:,jj) = W_H_inv(:,jj) ...
						./ (s - zF{jj}(ii));
					H_exp_inv(:,jj) = H_exp_inv(:,jj) ...
						.* (s - zF{jj}(ii));
				end
			end
		end

		ind_sc = ~isfinite(H_exp_inv);
		H_exp_inv(ind_sc) = 0;
		W_H_inv(ind_sc) = 0;
    end


%
%	PART II.	Main Process
%

	%
	%	Stacked Matrices for the System
	%
	H_exp = reshape(H_exp,n_freq*n_io,1);
	W_H	  = reshape(W_H  ,n_freq*n_io,1);

	if ind_Inv == 1,
		H_exp_inv = reshape(H_exp_inv,n_freq*n_io,1);
		W_H_inv	  = reshape(W_H_inv  ,n_freq*n_io,1);
	end

	sE = [];
	for ii=1:n_io,
		sE = [sE; ...
			s];
	end

	n_max = max(n_pole_R, max(n_zero_R) + 1);
	Phi = ones(n_freq,n_max);
	for ii = 2:n_max,
		Phi(:,ii) = s.^(ii-1);
	end

	%
	%	ARMA Parameters
	%
		%
		%	Phi_a
		%
		Phi_a = zeros(n_freq*n_io,n_pole_R);
		for ii = 1:n_io,
			Phi_a([1:n_freq]+(ii-1)*n_freq,:) = ...
				Phi(:,[n_pole_R:-1:1]);
		end

		%
		%	Phi_B
		%
		Phi_B = zeros(n_freq*n_io,sum(n_zero_R)+n_io);
		for ii=1:n_io,
			ind_B = (n_zero_R(ii)+1):-1:1;
			Phi_B([1:n_freq]+(ii-1)*n_freq, ...
				[1:(n_zero_R(ii)+1)] + ...
				ii-1+sum(n_zero_R(1:ii-1))) = ...
				Phi(:,ind_B);
		end

		if ind_SR == 1,
		%
		%	Phi_a_c & Phi_B_c
		%
		Phi_a_c = zeros(n_freq*n_io,n_pole_R);
		Phi_B_c = Phi_B;
		for ii=n_io:-1:1,
			n_c = size(TBa{ii},2);
			Phi_a_c([1:n_freq]+(ii-1)*n_freq, ...
				(n_pole_R+1-[n_c:-1:1])) = ...
				Phi(:,[n_c:-1:1]) * TBa{ii};
			Phi_B_c(:,n_zero_R(ii)+2+ ...
				sum(n_zero_R(1:ii-1))+(ii-1)- ...
				[n_c:-1:1]) = [];
		end
		end

	%
	%	Initial Parameter
	%
	sys_R = handles.sys_pool{handles.ind_sys};
	sys_F = sys_R;

	if ind_Fixed == 1,
		for jj=1:handles.n_i,
			for ii=1:handles.n_o,
				sys_F(ii,jj) = zpk(zF{ii+(jj-1)*n_o},pF, ...
					1);
				sys_R(ii,jj) = minreal(sys_R(ii,jj) ...
					/ sys_F(ii,jj));
            end
		end
	end
	[num,den] = tfdata(sys_R);

	theta_a = den{1,1}.';

	ind_srt = length(theta_a) - n_pole_R;
	if ind_srt > 1,
		if handles.ind_ck == 1,
	    disp('check coef of pole')
	    disp(theta_a(1:ind_srt-1));
	    end
	    theta_a = theta_a(ind_srt:end);
	end

    theta_B = [];
	for jj=1:n_i,
		for ii=1:n_o,
            theta_B_imsi = num{ii,jj}.' / theta_a(1);

            ind_srt = length(theta_B_imsi) ...
                - n_zero_R(ii+(jj-1)*n_o);
            if ind_srt > 1,
            	if handles.ind_ck == 1,
                disp('check coef of zero')
                disp(theta_B_imsi(1:ind_srt-1));
                end
                theta_B_imsi = theta_B_imsi(ind_srt:end);
            end

            if ind_SR == 0,
	            theta_B = [theta_B; ...
					theta_B_imsi];
			else,
				n_c = size(TBa{ii+(jj-1)*n_o},2);
	            theta_B = [theta_B; ...
					theta_B_imsi([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c])];
			end
		end
	end
	theta_a = theta_a / theta_a(1);
	theta_a(1) = [];

	theta = [theta_a;theta_B];

	%
	%	Initial Estimate
	%
	if ind_SR == 0,
		H_est = (Phi_B * theta_B) ...
			./ (sE.^n_pole_R + Phi_a * theta_a);
	else,
		H_est = (Phi_B_c * theta_B ...
			+ Phi_a_c * theta_a) ...
			./ (sE.^n_pole_R + Phi_a * theta_a);
	end

	eV = (H_exp-H_est).*W_H;
	if ind_Inv == 1,
		H_est_inv = 1 ./ H_est;
		eV = [eV;
			(H_exp_inv-H_est_inv).*W_H_inv];
	end
	ind_ck = ~isfinite(eV);
	ind_ck = logical(sign(sum(ind_ck.')));
	eV(ind_ck,:) = 0;

	eN = eV'*eV;
	eN_pool = [eN];

	eN_min = eN;

	if (ind_verb),
%	    clc
	    disp(['INITIAL ESTIMATE'])
		disp(['Current fit: ' num2str(eN)])
%		disp(['par-vector'])
%		disp(theta)
	end

	h_f_lm = figure;
	semilogy(0,eN_pool(1),'or')
	xlabel('Iteration #')
	ylabel('Error Norm')
	title('Levenberg-Marquardt Method')
	hold on

	%
	%	Levenberg-Marquardt
	%
		if ind_verb == 0,
		    options = optimset('Jacobian','on','Display','off');
		else,
		    options = optimset('Jacobian','on','Display','on');
		end

        if ind_Inv == 0,
            if ind_SR == 0,
        		thetaN = lsqnonlin(@mrtf_lm,theta,[],[],options,...
        		    H_exp, W_H, [], [], ...
    			n_io, n_pole_R, n_zero_R, ...
        		ind_Real, ind_SR, ind_Inv, ...
    			TBa, sE, ...
    			Phi_a, [], Phi_B, []);
            else,
        		thetaN = lsqnonlin(@mrtf_lm,theta,[],[],options,...
        		    H_exp, W_H, [], [], ...
    			n_io, n_pole_R, n_zero_R, ...
        		ind_Real, ind_SR, ind_Inv, ...
    			TBa, sE, ...
    			Phi_a, Phi_a_c, Phi_B, Phi_B_c);
            end
        else,
            if ind_SR == 0,
        		thetaN = lsqnonlin(@mrtf_lm,theta,[],[],options,...
        		    H_exp, W_H, H_exp_inv, W_H_inv, ...
        			n_io, n_pole_R, n_zero_R, ...
            		ind_Real, ind_SR, ind_Inv, ...
        			TBa, sE, ...
        			Phi_a, [], Phi_B, []);
            else,
        		thetaN = lsqnonlin(@mrtf_lm,theta,[],[],options,...
        		    H_exp, W_H, H_exp_inv, W_H_inv, ...
        			n_io, n_pole_R, n_zero_R, ...
            		ind_Real, ind_SR, ind_Inv, ...
        			TBa, sE, ...
        			Phi_a, Phi_a_c, Phi_B, Phi_B_c);
            end
        end

    figure(h_f_lm);
	semilogy(eN_pool(2:end),'o')
	semilogy(size(eN_pool,1)+1,eN_pool(end),'or')
	hold off

%
%	PART III.	Post Process
%
	%
	%	Parameter Vectors
	%
	theta = thetaN;

	theta_a_R = [1;theta(1:n_pole_R)];
   	theta_a_R = apolystab(theta_a_R,ind_Real);
    theta_B_R = cell(n_io,1);
	if ind_SR == 0,
		for jj=1:n_i,
			for ii=1:n_o,
				theta_B_R{ii+(jj-1)*n_o} = ...
					theta([1:(n_zero_R(ii+(jj-1)*n_o)+1)] + ...
						sum(n_zero_R(1:(ii-1+(jj-1)*n_o))) + ...
						(ii-1+(jj-1)*n_o) + ...
						n_pole_R);
			end
		end
	else,
		ind_th = n_pole_R;
		for jj=1:n_i,
			for ii=1:n_o,
				n_c = size(TBa{ii+(jj-1)*n_o},2);
				theta_B_R{ii+(jj-1)*n_o} = ...
					[theta([1:n_zero_R(ii+(jj-1)*n_o)+1-n_c] + ...
						ind_th);
						TBa{ii+(jj-1)*n_o}*...
							theta_a_R(n_pole_R+2-[n_c:-1:1]);
					];
				ind_th = ind_th + ...
					n_zero_R(ii+(jj-1)*n_o)+1-n_c;
			end
		end
	end

	%
	%	Rational Polynomial Model
	%
	for jj=1:n_i,
		for ii=1:n_o,
			num{ii,jj} = theta_B_R{ii+(jj-1)*n_o}.';
			den{ii,jj} = theta_a_R.';
		end
	end
	sys_R = tf(num,den);

	if ind_Fixed == 1,
		[z,p,k] = zpkdata(sys_R);
		for jj=1:n_i,
			for ii=1:n_o,
        		p{ii,jj} = [p{ii,jj}; ...
                    handles.fixed_pole];
				z{ii,jj} = [z{ii,jj}; ...
                        handles.fixed_zero{ii+(jj-1)*n_o}];
			end
		end
		sys_R = zpk(z,p,k);
	else,
		handles.fixed_pole = [];
		for ii = 1:handles.n_io,
			handles.fixed_zero{ii} = [];
		end
	end

	%
	%	Data Base & Graphic
	%
	handles.ind_sys = length(handles.sys_pool) + 1;
    handles.sys_pool{handles.ind_sys} = ...
 		sys_R;

	str_ind_sys = [ ...
		num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);
	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	    pz_window(handles);
	end

		%   Update handles structure
		guidata(hObject, handles);


% --- Executes on button press in pushbutton_ReC.
function pushbutton_ReC_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ReC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton_ReC

%	warning off
%	echo off

%   Basic Parameter Setting
	im = sqrt(-1);

	freq  = handles.freq(handles.ind_freq);
    omega = freq*2*pi;
    s	  = im*omega;

    nm	 = handles.n_o;
	nu	 = handles.n_i;
    ns   = handles.n_io;
    nf   = length(freq);

%   Experiment Data
    Hyf_m = handles.H_exp(handles.ind_freq,:);
	sig_H = handles.W_H(handles.ind_freq,:).^(1/2);

%   Problem Type
%	index_mp = 0;
%	if get(handles.togglebutton_Inv,'Value') > 0,
%	index_mp = 1;
%	end
%	index_sr = handles.ind_SR;


%
%	PART II.	Standard Form
%
	sys_o = handles.sys_pool{handles.ind_sys};
	sys_om = minreal(sys_o);

%   Model Data
	[b,a] = tfdata(sys_om);

%	Residue Model
	for ii=1:nm,
	for jj=1:nu,
		bij = b{ii,jj};
		aij = a{ii,jj};
	if 0,
	for kk = 1:length(bij),
		fprintf(1,'% 10.5e & % 10.5e :',real(bij(kk)),imag(bij(kk)))
		fprintf(1,'% 10.5e & % 10.5e :',real(aij(kk)),imag(aij(kk)))
		fprintf(1,'\n')
	end
	end
	end
	end


%
%	PART III.	Modal Parameter Extraction
%
	lambda = [];
	gamma  = [];
	r_p    = [];
	k_p    = [];

	for ii=1:nm,
	for jj=1:nu,
		a{ii,jj} = apolystab(a{ii,jj},1);
		[r,p,k]  = residue(b{ii,jj},a{ii,jj});

		wD	 = abs(imag(p));
		wN   = abs(p);
		fN   = wN / (2*pi);
		zeta =-real(p)./wN;

		lambda = [lambda, ...
			p];
		gamma_p = -real(p) - im*wN.^2.*(1-2*zeta.^2)./(2*imag(p));
		gamma = [gamma, ...
			gamma_p];
		r_p = [r_p, ...
			r./gamma_p];
		k_p = [k_p, ...
			k];
	end
	end

%	phi & psi & mu
	n_mode = ceil(length(wD)/2);
	r_p = reshape(r_p,ns*2,n_mode);

	%	Optimization
	phi = zeros(nm,n_mode);
	psi = zeros(1,n_mode);
	mui = zeros(1,n_mode);

	theta0 = zeros(ns,n_mode);
	theta1 = zeros(ns,n_mode);

%		ii = n_mode;							% mode #; 1:n_mode
	for ii = n_mode:-1:1,
		ind_off = (ii-1)*2;
		y_cm = r_p(:,ii);

		%
%		y_cm = log(y_cm);
%		phi_cm01 = zeros(nm*nu*2,nm);
%		phi_cm02 = zeros(nm*nu*2,nu*2);
%		for jj=1:nu*2,
%			phi_cm01([1:nm]+(jj-1)*nm,:) = eye(nm);
%			phi_cm02([1:nm]+(jj-1)*nm,jj) = im;
%		end
%
%		phi_cm = [phi_cm01,phi_cm02];
%		theta_cm = inv(phi_cm.'*phi_cm)*phi_cm.'*y_cm;

		phi0 = y_cm;
		phi0([1:2:ns*2]) = phi0([1:2:ns*2]) ./ phi0([2:2:ns*2]);
		phi0([2:2:ns*2]) = phi0([2:2:ns*2]) ./ phi0([2:2:ns*2]);

		phi0 = reshape(phi0,2,ns);
%		w_d = log(phi1 ./ (phi1(:,1)*ones(1,ns))) * 180/pi;

		phi0 = mean(phi0.').';
		phi0 = real(phi0);
		phi(:,ii) = phi0 / sqrt(sum(phi0.^2));

		psi(1,ii) = ...
		mean([ ...
		y_cm([1:nm]   ,1) ./ phi(:,ii) - phi(1,ii) - phi(2,ii)
		y_cm([1:nm]+ns,1) ./ phi(:,ii) - phi(1,ii) - phi(2,ii) ...
		]);
		psi = real(psi);

		mui(1,ii) = ...
		mean([ ...
		-(phi(:,ii) ./ y_cm([1:nm]+nm   ,1)) / phi(1,ii)
		-(phi(:,ii) ./ y_cm([1:nm]+nm+ns,1)) / phi(1,ii) ...
		]);

		mui = abs(real(mui));

		theta0(:,ii) = [phi(:,ii);psi(ii);mui(ii)];
		eN_pool = [];
		r2nm00(theta0(:,ii),y_cm);

		lb = [-Inf, -Inf, -Inf,  0];
		ub = [ Inf,  Inf,  Inf,  Inf];

		options = optimset('Jacobian','on','Display','on');
%		options = optimset('Jacobian','on','Display','on','TolX',1e-15);
%		options = optimset('Jacobian','on','Display','on', ...
%			'MaxFunEvals',1e4,'MaxIter',1e6,'TolX',1e-15);
%		options = optimset('Jacobian','on','Display','on', ...
%			'TolX',1e-15,'TolFun',1e-18);
%		options = optimset('Jacobian','on','Display','on', ...
%			'MaxFunEvals',1e4,'TolFun',1e-9,'MaxIter',1e6);
		theta1(:,ii) = lsqnonlin(@r2nm00,theta0(:,ii),lb,ub,options,...
			y_cm);
		e_tr{ii} = eN_pool;
		pause

	end

	figure(1)
	for ii=1:n_mode,
		subplot(2,2,ii)
		eN_pool = e_tr{ii};
		[min_e, ind_e] = min(eN_pool(:,1));
		theta1m(:,ii) = eN_pool(ind_e,[2:5]).';

%		plot(eN_pool(:,1))
		semilogy(eN_pool(:,1))
		hold on
		semilogy(1,eN_pool(1,1),'ro')
		semilogy(size(eN_pool,1),eN_pool(end,1),'ro')
		hold off
		xlabel([num2str(min_e),' @ ',num2str(ind_e)])
		ylabel([num2str(ii),'th Mode'])
	end


	%	Optimization w/ Constraint
	eN_pool = [];
%	r2nm10(theta0(:),r_p);
	r2nm10(theta1m(:),r_p);
%	r2nm10(theta1(:),r_p);

	lb = [-Inf, -Inf, -Inf,  0];	lb = [lb,lb,lb,lb];
	ub = [ Inf,  Inf,  Inf,  Inf];	ub = [ub,ub,ub,ub];

	options = optimset('Jacobian','on','Display','on');
%	options = optimset('Jacobian','on','Display','on', ...
%		'MaxFunEvals',1e6);
%	options = optimset('Jacobian','on','Display','on', ...
%		'MaxFunEvals',1e4,'MaxIter',1e6);
%		'MaxFunEvals',1e4,'MaxIter',1e6,'TolX',1e-15);
	theta2 = fmincon(@r2nm10,theta1m(:),[],[],[],[],lb,ub, ...
		@r2nm10c,options, ...
		r_p);

	figure(2)
	[min_e, ind_e] = min(eN_pool(100:end,1));	ind_e = ind_e + 100;
	semilogy(eN_pool(:,1))
	hold on
	semilogy(1,eN_pool(1,1),'ro')
	semilogy(ind_e,min_e,'ro')
	semilogy(size(eN_pool,1),eN_pool(end,1),'ro')
	hold off
	xlabel([num2str(eN_pool(1,1)), ', ', ...
		num2str(min_e), ', ', ...
		num2str(eN_pool(end,1)), ...
		' @ ',num2str(ind_e)])
	ylabel([num2str(ii),'th Mode'])

	theta2 = reshape(theta2,4,4);
	phi2 = theta2(1:3,:);

	phi2'*phi2
	theta2

	if 0,
	theta2m = reshape(eN_pool(ind_e,2:end),4,4);
	phi2m = theta2m(1:3,:);

	phi2m'*phi2m
	theta2m
	end


%
%	PART III.	Normal Modal T.F. Model
%
	n_th = (4+nm);
	theta30 = zeros(n_th*n_mode,1);
	theta30(1:n_mode) = (wN(1:2:end) + wN(2:2:end))/2;
	theta30([1:n_mode]+n_mode) = (zeta(1:2:end) + zeta(2:2:end))/2;
	theta30([1:n_mode]+n_mode*2) = theta2(nm+2,:).';
	theta30([1:n_mode*(nm+1)]+n_mode*3) = phi2(:);

	%	Optimization w/o Constraint
	eN_pool = [];
	r2nm20(theta30, ...
		Hyf_m,sig_H, ...
		n_mode, nm, ...
		s,nf);

	lb = [zeros(1,n_mode), zeros(1,n_mode), zeros(1,n_mode), ...
		-Inf*ones(1,n_mode), -Inf*ones(1,n_mode), -Inf*ones(1,n_mode)];
	ub = [Inf*ones(1,n_mode), Inf*ones(1,n_mode), Inf*ones(1,n_mode), ...
		Inf*ones(1,n_mode), Inf*ones(1,n_mode), Inf*ones(1,n_mode)];

%	options = optimset('Jacobian','on','Display','on');
	options = optimset('Jacobian','on','Display','on', ...
		'TolX',1e-15);
%		'MaxFunEvals',1e6);
%	options = optimset('Jacobian','on','Display','on', ...
%		'MaxFunEvals',1e4,'MaxIter',1e6);
%		'MaxFunEvals',1e4,'MaxIter',1e6,'TolX',1e-15);
	theta31 = lsqnonlin(@r2nm20,theta30,lb,ub,options,...
		Hyf_m,sig_H, ...
		n_mode, nm, ...
		s,nf);

	figure(3)
	[min_e, ind_e] = min(eN_pool(100:end,1));	ind_e = ind_e + 100;
	semilogy(eN_pool(:,1))
	hold on
	semilogy(1,eN_pool(1,1),'ro')
	semilogy(ind_e,min_e,'ro')
	semilogy(size(eN_pool,1),eN_pool(end,1),'ro')
	hold off
	xlabel([num2str(eN_pool(1,1)), ', ', ...
		num2str(min_e), ', ', ...
		num2str(eN_pool(end,1)), ...
		' @ ',num2str(ind_e)])
	ylabel([num2str(ii),'th Mode'])

	%
	wN31 = theta31(1:n_mode);
	zeta31 = theta31([1:n_mode]+n_mode);
	mu31 = theta31([1:n_mode]+n_mode*2);
	phi31 = theta31([1:n_mode*(nm+1)]+n_mode*3);
	phi31 = reshape(phi31,(nm+1),n_mode);

	phi31'*phi31


	%	Optimization w/ Constraint
%	warning off
%	echo off
	eN_pool = [];
	r2nm21(theta31, ...
		Hyf_m,sig_H, ...
		n_mode, nm, ...
		s,nf);

	options = optimset('Jacobian','on','Display','on');
%	options = optimset('Jacobian','on','Display','on', ...
%		'MaxFunEvals',1e6);
%	options = optimset('Jacobian','on','Display','on', ...
%		'MaxFunEvals',1e4,'MaxIter',1e6);
%		'MaxFunEvals',1e4,'MaxIter',1e6,'TolX',1e-15);
	theta32 = fmincon(@r2nm21,theta31,[],[],[],[],lb,ub, ...
		@r2nm21c,options, ...
		Hyf_m,sig_H, ...
		n_mode, nm, ...
		s,nf);


	figure(5)
	[min_e, ind_e] = min(eN_pool(100:end,1));	ind_e = ind_e + 100;
	semilogy(eN_pool(:,1))
	hold on
	semilogy(1,eN_pool(1,1),'ro')
	semilogy(ind_e,min_e,'ro')
	semilogy(size(eN_pool,1),eN_pool(end,1),'ro')
	hold off
	xlabel([num2str(eN_pool(1,1)), ', ', ...
		num2str(min_e), ', ', ...
		num2str(eN_pool(end,1)), ...
		' @ ',num2str(ind_e)])
	ylabel([num2str(ii),'th Mode'])

	wN32 = theta32(1:n_mode);
	zeta32 = theta32([1:n_mode]+n_mode);
	mu32 = theta32([1:n_mode]+n_mode*2);
	phi32 = theta32([1:n_mode*(nm+1)]+n_mode*3);
	phi32 = reshape(phi32,(nm+1),n_mode);

	phi32'*phi32


%
%	PART IV.	Re Construction
%
	wN40   = wN31;
	zeta40 = zeta31;
	mu40   = mu31;
	phi40  = phi31;

	lambda40 = -zeta40.*wN40 + im*wN40 .* sqrt(1-zeta40.^2);
	alpha40  = - im./(2*wN40.*sqrt(1-zeta40.^2));
	beta40   = 2*sqrt(1-zeta40.^2)./zeta40;
	gamma40 = zeta40.*wN40 - im*wN40.*(1-2*zeta40.^2) ./ sqrt(1-zeta40.^2) / 2;

	Anm = zeros(2*n_mode);
	Bnm = zeros(2*n_mode,nu);
	Cnm = zeros(nm*3,2*n_mode);
	Dnm = zeros(nm*3,nu);

	for ii=1:n_mode,
		Anm([1:2]+(ii-1)*2,[1:2]+(ii-1)*2) = ...
			[real(lambda40(ii)), imag(lambda40(ii));
			-imag(lambda40(ii)), real(lambda40(ii))];
		Bnm(1+(ii-1)*2,[1:2]) = ...
			[-sum(phi40(:,ii)),phi40(1,ii)/mu40(ii)]*2;
		Cnm([1:nm*3],[1:2]+(ii-1)*2) = ...
			[phi40(1:nm,ii)*[real(alpha40(ii)),imag(alpha40(ii))];
			 phi40(1:nm,ii)*[real(beta40(ii)),imag(beta40(ii))];
			-phi40(1:nm,ii)*[real(gamma40(ii)),imag(gamma40(ii))];
			];
		Dnm([1:nm]+2*nm,nu) = Dnm([1:nm]+2*nm,nu) + ...
			phi40(1:nm,ii)/mu40(ii)*phi40(1,ii);
	end

	sys_nm = ss(Anm,Bnm,Cnm,Dnm);

	% Update handles structure
	guidata(hObject, handles);


	global eN_pool;

%	warning off
%	echo off

	%
	%	Measurement Data Set
	%
		nm	 = handles.n_o;
		nu	 = handles.n_i;
	    ns   = handles.n_io;
	    nf   = length(freq);

		freq = handles.freq(handles.ind_freq)*2*pi;

		Hyf  = handles.H_exp(handles.ind_freq,:);
		W    = handles.W_H(handles.ind_freq,:);

		Hyf_inv = 1./Hyf;
		Hyf_inv(~isfinite(Hyf_inv)) = 0;
		W_inv 	= handles.W_H_inv(handles.ind_freq,:);

    %
    %	Initial Condition
	%
		[sys_am,str] = model_am;

		m1 = str.Ms(1,1);
		m2 = str.Ms(2,2);
		c2 = str.Cs(2,2);
		c1 = str.Cs(1,1)+str.Cs(1,2);
%		c1 = str.Cs(1,1);
		k2 = str.Ms(2,2);
		k1 = str.Ks(1,1)+str.Ks(1,2);
%		k1 = str.Ks(1,1);
		sen_fd = 1/110e-3;
%		sen_fd = 1e10;
		ga = 9.80665;

%		th0 = [m1,m2,c1,c2,k1,k2,sen_fd,ga];
%		th0 = [m1,m2,c1,c2,k1,k2,sen_fd];
%		th0 = [m1,c1,c2,k1,k2,sen_fd];
%		th0 = [m1,c1,c2,k1,k2];
%		th0 = [c1,c2,k1,k2];
        th0 = [k1,k2];

	%
	%	Levenberg-Marquardt
	%
	    options = optimset('Display','off');
%	    options = optimset('Jacobian','on','Display','off');

        lb = 0*ones(size(th0));
        ub = Inf*ones(size(th0));

		eV = mrtf_str1(th0, ...
			freq,Hyf,W, Hyf_inv, W_inv, ...
            m1,m2,c1,c2,ga,sen_fd);
%        m1,m2,ga,sen_fd);
		th = lsqnonlin(@mrtf_str1,th0,lb,ub,options,...
			freq,Hyf,W, Hyf_inv, W_inv, ...
            m1,m2,c1,c2,ga,sen_fd);
%        m1,m2,ga,sen_fd);

        figure(5)
		semilogy(0,eN_pool(1),'or')
		xlabel('Iteration #')
		ylabel('Error Norm')
		hold on
		semilogy(eN_pool(2:end),'o')
		semilogy(size(eN_pool,1),eN_pool(end),'or')
		hold off

        ind_th = 4;
%	m1		= th(1);
%	m2      = th(2);
%	c1      = th(3-ind_th);
%	c2      = th(4-ind_th);
	k1      = th(5-ind_th);
	k2      = th(6-ind_th);
%	sen_fd  = th(7-ind_th);
%	ga      = th(8);

	%
	%	Structure Model
	%

	MsP = [m1,m2];
	Ms = diag(MsP);						iMs = inv(Ms);
%	Cs = [c1, -c2;
	Cs = [c1+c2, -c2;
		 -c2, c2];
%	Ks = [k1, -k2;
	Ks = [k1+k2, -k2;
		 -k2, k2];

%	Bu = [-m1/ga,sen_fd;
%		  -m2/ga,0];

	%
	%	State Space Model
	%

	Ao = [zeros(2),  eye(2);
		 -iMs*Ks,  -iMs*Cs];
	Bo = [zeros(2);
		 -[1;1]*ga, [sen_fd/m1;0]];
	Co = Ao([3:4],:)/ga;
	Do = Bo([3:4],:)/ga;
    Do([1:2],1) = 0;

	sys_o = ss(Ao,Bo,Co,Do);
%	sys_o = sys_am;


%
%	FINE
%


% --- Executes on button press in togglebutton_SR.
function togglebutton_SR_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_SR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_SR

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
    set(hObject,'FontWeight','bold');
    set(hObject,'ForegroundColor',[1,1,0]);
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
    set(hObject,'FontWeight','normal');
    set(hObject,'ForegroundColor',[0,0,0]);
end


% --- Executes on button press in pushbutton_MR.
function pushbutton_MR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	%
	%	Basic Parameter Setting
	%
	im = sqrt(-1);

	n_o      = handles.n_o;
	n_i      = handles.n_i;
	n_io     = handles.n_io;

	n_pole	 = handles.n_pole;
	n_zero	 = handles.n_zero;
	n_pole_F = 0;
	n_zero_F = zeros(n_io,1);
    n_pole_R = n_pole - n_pole_F;
    n_zero_R = n_zero - n_zero_F;

	pF		 = handles.fixed_pole;
	zF		 = handles.fixed_zero;

	%
	%	Option : Real Flag
	%
	ind_Real = 0;
    if get(handles.togglebutton_Real,'Value') > 0,
	    ind_Real = 1;
	end

	%
	%	Option : Fixed
	%
	ind_Fixed = 0;
    if get(handles.togglebutton_Fixed_PZ,'Value') > 0,
	    ind_Fixed = 1;

	    n_pole_F = length(pF);
		for ii=1:n_io,
			n_zero_F(ii) = length(zF{ii});
	    end
	    n_pole_R = n_pole - n_pole_F;
	    n_zero_R = n_zero - n_zero_F;
    end

	%
	%	Estimated System
	%
	sys_R = handles.sys_pool{handles.ind_sys};
	sys_F = sys_R;

	if ind_Fixed == 1,
		for jj=1:handles.n_i,
			for ii=1:handles.n_o,
				sys_F(ii,jj) = zpk(zF{ii+(jj-1)*n_o},pF, ...
					1);
				sys_R(ii,jj) = minreal(sys_R(ii,jj) ...
					/ sys_F(ii,jj));
            end
		end
	end
	[num,den] = tfdata(sys_R);

	%
	%	Minimum System Realization
	%
	theta_a = den{1,1}.';
   	theta_a = apolystab(theta_a,ind_Real);

	ind_srt = length(theta_a) - n_pole_R;
	if ind_srt > 1,
		if handles.ind_ck == 1,
		disp('check coef of pole')
		disp(theta_a(1:ind_srt-1));
		end
		theta_a = theta_a(ind_srt:end);
	end

	for jj=1:n_i,
		for ii=1:n_o,
			theta_B_ij = num{ii,jj}.' / theta_a(1);

            ind_srt = length(theta_a) - n_pole_R;

            if ind_srt > 1,
            	if handles.ind_ck == 1,
                disp('check coef of pole')
                disp(theta_a(1:ind_srt-1));
                end
                theta_a = theta_a(ind_srt:end);
            end

			theta_B_ij = apolystab(theta_B_ij,ind_Real);
            num{ii,jj} = theta_B_ij.';
		end
	end
	theta_a = theta_a / theta_a(1);

	%
	%	Rational Polynomial Model
	%
	for jj=1:n_i,
		for ii=1:n_o,
			den{ii,jj} = theta_a.';
		end
	end
	sys_R = tf(num,den);

	if ind_Fixed == 1,
		[z,p,k] = zpkdata(sys_R);
		for jj=1:n_i,
			for ii=1:n_o,
        		p{ii,jj} = [p{ii,jj}; ...
                    handles.fixed_pole];
				z{ii,jj} = [z{ii,jj}; ...
                        handles.fixed_zero{ii+(jj-1)*n_o}];
			end
		end
		sys_R = zpk(z,p,k);
	else,
		handles.fixed_pole = [];
		for ii = 1:handles.n_io,
			handles.fixed_zero{ii} = [];
		end
	end

	%   Scaling
    ind_f_eval = sum(handles.freq <= handles.f_eval);

	[z,p,k] = zpkdata(sys_R);
    for ii=1:handles.n_io,
	    sys_est= zpk(z{ii},p{ii},1);
	    H_est  = freqresp(sys_est,handles.freq(ind_f_eval)*2*pi);
	    H_exp  = handles.H_exp(ind_f_eval,ii);

%	    k(ii) = abs(h_ref)/abs(h_imsi);
	    k(ii) = H_exp/H_est;
    end
	sys_R = zpk(z,p,k);

	%
	%	Data Base & Graphic
	%
    handles.ind_sys = handles.ind_sys + 1;
    handles.sys_pool{handles.ind_sys} = ...
 		sys_R;

	str_ind_sys = [ ...
		num2str(handles.ind_sys), ...
		'/', ...
		num2str(length(handles.sys_pool))];
	set(handles.text_ind_sys, 'String', str_ind_sys);
	n_error_norm(handles);

	draw_n(handles);
	draw_frf(handles);
    try,
	    pz_window(handles);
	end

		%   Update handles structure
		guidata(hObject, handles);


% --- Executes on button press in pushbutton_Cont_In.
function pushbutton_Cont_In_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cont_In (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if 0,
	freq_range_0 = handles.freq(diff(handles.ind_freq)==1);
	freq_range_1 = handles.freq(diff(handles.ind_freq)==-1);
end

	data_input = datain;

if 0,
    handles.freq = data_input.freq;
    handles.H_exp = data_input.H_exp;
    handles.W_H = data_input.W_H;
    handles.W_H_inv = data_input.W_H_inv;

    handles.n_o = data_input.n_o;
    handles.n_i = data_input.n_i;
    handles.n_io = data_input.n_i*data_input.n_o;
end

    handles.n_pole = data_input.n_pole;
    handles.n_zero = data_input.n_zero*ones(handles.n_io,1);

    handles.ind_SR_o = data_input.ind_SR_o;
    handles.ind_SR_i = data_input.ind_SR_i;

if 0,
	handles = data_pre(handles,1);

	iind1 = 0;
	if freq_range_1(1) < freq_range_0(1),
		handles.ind_freq = freq <= freq_range_1(1);
		iind1 = 1;
	end
	for ii=1:(length(freq_range_1)-iind1),
		handles.ind_freq = handles.ind_freq | ...
			((freq > freq_range_0(ii)) && ...
			 (freq <= freq_range_1(ii)));
	end
	if freq_range_1(end) < freq_range_0(end),
		handles.ind_freq = handles.ind_freq | ...
			(freq > freq_range_0(end));
	end
end

	set(handles.edit_n_Poles,'String',...
    	{num2str(handles.n_pole)});
	set(handles.edit_n_Zeros,'String',...
    	{num2str(handles.n_zero(handles.ind_f))});

    guidata(hObject, handles);


% --- Executes on button press in togglebutton_Fixed_PZ.
function togglebutton_Fixed_PZ_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_Fixed_PZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_Fixed_PZ

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
    set(hObject,'FontWeight','bold');
    set(hObject,'ForegroundColor',[1,1,0]);
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
    set(hObject,'FontWeight','normal');
    set(hObject,'ForegroundColor',[0,0,0]);
end


% --- Executes on button press in pushbutton_MP.
function pushbutton_MP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_Real.
function togglebutton_Real_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_Real (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_Real

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
    set(hObject,'FontWeight','bold');
    set(hObject,'ForegroundColor',[1,1,0]);
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
    set(hObject,'FontWeight','normal');
    set(hObject,'ForegroundColor',[0,0,0]);
end


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file = uigetfile('*.mat');
if ~isequal(file, 0)
%	pos_size = get(handles.fig_mfdid,'Position');
%	lb(file);
	data_input = datain(file);

    handles.freq = data_input.freq;
    handles.H_exp = data_input.H_exp;
    handles.W_H = data_input.W_H;
    handles.W_H_inv = data_input.W_H_inv;

    handles.n_o = data_input.n_o;
    handles.n_i = data_input.n_i;
    handles.n_io = data_input.n_i*data_input.n_o;

    handles.n_pole = data_input.n_pole;
    handles.n_zero = data_input.n_zero*ones(handles.n_io,1);

    handles.ind_SR_o = data_input.ind_SR_o;
    handles.ind_SR_i = data_input.ind_SR_i;

    handles = data_pre(handles,0);

    guidata(hObject, handles);

end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.fig_mfdid)


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selection = questdlg(['Close ' get(handles.fig_mfdid,'Name') '?'],...
                     ['Close ' get(handles.fig_mfdid,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.fig_mfdid)


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current position of the GUI from the handles structure
% to pass to the modal dialog.
pos_size = get(handles.fig_mfdid,'Position');
% Call modaldlg with the argument 'Position'.
user_response = AboutMFDID;


function load_listbox(dir_path, handles)
	cd (dir_path)
	dir_struct = dir(dir_path);
	[sorted_names,sorted_index] = sortrows({dir_struct.name}');
	handles.file_names = sorted_names;
	handles.is_dir = [dir_struct.isdir];
	handles.sorted_index = [sorted_index];
	guidata(handles.fig_mfdid,handles)
	set(handles.popupmenu_MDM_L,'String',handles.file_names,...
	    'Value',1)
% 	set(handles.text1,'String',pwd)
