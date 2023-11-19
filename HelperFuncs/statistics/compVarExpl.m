function g = compVarExpl(filename)

videoFile               = [filename '.mpg'];
procFile                = [filename '_proc.mat'];

if ~exist(procFile,'file')
    warning('No processed video data file found ')
end

load(procFile)

%% Initialize video reading:
cap                     = py.cv2.VideoCapture(videoFile);

fps                     = cap.get(int8(5)); %#https://docs.opencv.org/2.4/modules/highgui/doc/reading_and_writing_images_and_video.html#videocapture-get
nframes                 = cap.get(int8(7));
duration                = nframes*1./fps;
widthPixels             = cap.get(int8(3));
heightPixels            = cap.get(int8(4));

%% try to read some frames until it works, then reset: 
while cap.isOpened
    ret     = cap.read;
    mt      = cell(ret(1));
    if mt{1}
        iFr=1;
        break;
    end
end
cap.set(int8(1),0);

%%

h.nX = {[720]};
h.nY = {[576]};
h.ROI = {{[1 1 179 143]}};
ns = 4;

subsamplepix = 143*179;

nselec = 1000;
selecframes = round(nframes/2);
selecframes = selecframes:selecframes+1000-1;

moviedata       = NaN(nselec,subsamplepix);

cap.set(int8(1),selecframes(1)-1);

for iF = 1:nselec
    iF
    try
        ret     = cap.read;
        mt      = cell(ret(2));
        x       = mt{1};
        data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
        data    = data(1:3:end);
        im      = reshape(data,[h.nX{1} h.nY{1}])';
        
        im      = squeeze(mean(mean(reshape(single(im(1:floor(h.nY{1}/ns)*ns,1:floor(h.nX{1}/ns)*ns,1)),...
            ns, floor(h.nY{1}/ns), ns, floor(h.nX{1}/ns)), 1),3));
        im = im(proc.wpix{1});
        
%         ims{1}(np(k) + [1:npix(k)], t) = im(h.wpix{k}(:));

        moviedata(iF,:) = im(:);
        
    catch
        
    end
    
end

%%

figure;
% imagesc(flipud(reshape(moviedata(iF,:),143,179)'),[0 256])
subplot 151 
exframemot = abs(diff(moviedata([1 2],:),[],1));
exframemot = flipud(reshape(exframemot,143,179)');
imagesc(exframemot,[-5 5])
title(nanvar(exframemot(:)))

subplot 152
avgmot          = reshape(proc.avgmotion,144,180);
avgmot          = avgmot(proc.wpix{1});
avgmot          = flipud(reshape(avgmot,143,179)');
imagesc(avgmot,[-5 5])
title(nanvar(avgmot(:)))

subplot 153
subtr1 = exframemot-avgmot;
imagesc(subtr1,[-5 5])
title(nanvar(subtr1(:)))

subplot 154
iComp = 500;
temp            = proc.motSVD{1}(selecframes(1)+1,1:iComp) * proc.uMotMask{1}(:,1:iComp)';
temp            = flipud(reshape(temp,143,179)');
imagesc(temp,[-5 5])
title(nanvar(temp(:)))

subplot 155
subtr2 = exframemot-avgmot-temp;
imagesc(subtr2,[-5 5])
title(nanvar(subtr2(:)))

%% subtract mean motion:
moviemotdata    = [NaN(1,subsamplepix); abs(diff(moviedata,[],1))];
% imot = bsxfun(@minus, abs(diff(ims{z},1,2)), avgmot{z});

avgmot          = reshape(proc.avgmotion,144,180);
avgmot          = avgmot(proc.wpix{1});
moviemotdata    = moviemotdata - repmat(avgmot(:)',nselec,1);

%%
for iComp = 1:500
    iComp
    temp            = proc.motSVD{1}(selecframes,1:iComp) * proc.uMotMask{1}(:,1:iComp)'; 
%     temp = temp+repmat(proc.avgmotion
%     g(iComp)        = nanvar(temp(:)) / nanvar(moviemotdata(:));
    g(iComp)        = 1 - nanvar(moviemotdata(:) - temp(:)) / nanvar(moviemotdata(:));
end

end