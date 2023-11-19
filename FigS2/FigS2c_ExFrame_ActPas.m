%% Script that shows example frames from active and passive task versions:
% MOL (C) 2023
% as reported in Oude Lohuis et al. 2023 Nat Neurosci
% "Triple dissociation of auditory, visual and motor processing in primary visual cortex"

%# DATA NOT INCLUDED IN THE REPOSITORY! 

startover

%% 
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\3ActivePassive';

params.filename_pas = 'F:\Data\CHDET\RawData\2045\2021-04-28_19-13-59\ChangeDetectionConflict\VT1.0001.mpg';
params.filename_act = 'F:\Data\CHDET\RawData\2045\2021-04-28_19-13-59\ChangeDetectionConflict2\VT1.0002.mpg';

%%
widthPixels = 720;
heightPixels = 576;

%% Figure of video frames active/passive:

params.exFr_act = 1149;

params.exFr_pas = 1000;

figure; subplot(1,2,1)
cap                     = py.cv2.VideoCapture(params.filename_act);
cap.set(int8(1), params.exFr_act-1);
ret     = cap.read();
mt      = cell(ret(2));
x       = mt{1};
data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
data    = data(1:3:end);
data    = reshape(data,[widthPixels heightPixels]);
imshow(uint8(data));
set(gca,'YDir','normal')
title('','FontSize',20)

subplot(1,2,2)
cap                     = py.cv2.VideoCapture(params.filename_pas);
cap.set(int8(1), params.exFr_pas-1);
ret     = cap.read();
mt      = cell(ret(2));
x       = mt{1};
data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
data    = data(1:3:end);
data    = reshape(data,[widthPixels heightPixels]);
imshow(uint8(data));
set(gca,'YDir','normal')
title('','FontSize',20)

filename = sprintf('ExFrame_ActivePassive.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%% 
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\3ActivePassive';

params.filename_pas = 'F:\Data\CHDET\RawData\2045\2021-04-29_15-55-10\ChangeDetectionConflict\VT1.0001.mpg';
params.filename_act = 'F:\Data\CHDET\RawData\2045\2021-04-29_15-55-10\ChangeDetectionConflict2\VT1.0002.mpg';

%% Figure of video frames active/passive:

params.exFr_act = 10847;

params.exFr_pas = 31040;

figure; subplot(1,2,1)
cap                     = py.cv2.VideoCapture(params.filename_act);

%try to read some frames until it works, then reset: 
while cap.isOpened
    ret     = cap.read;
    mt      = cell(ret(1));
    if mt{1}
        break;
    end
end
% cap.set(int8(1),0);

cap.set(int8(1), params.exFr_act-1);

ret     = cap.read();
mt      = cell(ret(2));
x       = mt{1};
data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
data    = data(1:3:end);
data    = reshape(data,[widthPixels heightPixels]);
imshow(uint8(data));
set(gca,'YDir','normal')
title('','FontSize',20)

subplot(1,2,2)
cap                     = py.cv2.VideoCapture(params.filename_pas);

%try to read some frames until it works, then reset: 
while cap.isOpened
    ret     = cap.read;
    mt      = cell(ret(1));
    if mt{1}
        break;
    end
end
% cap.set(int8(1),0);

cap.set(int8(1), params.exFr_pas-1);
ret     = cap.read();
mt      = cell(ret(2));
x       = mt{1};
data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
data    = data(1:3:end);
data    = reshape(data,[widthPixels heightPixels]);
imshow(uint8(data));
set(gca,'YDir','normal')
title('','FontSize',20)

filename = sprintf('ExFrame_ActivePassive_v2.eps');
export_fig(fullfile(params.savedir,filename),gcf);

