% Read in a standard MATLAB color demo image.
folder = fullfile(matlabroot, '\toolbox\images\imdemos');
baseFileName = 'peppers.png';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
rgbImage = imread(fullFileName);
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows columns numberOfColorBands] = size(rgbImage);
% Display the original color image.
subplot(2, 2, 1);
imshow(rgbImage, []);
axis on;
title('Original Color Image');
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
smallSubImage = imcrop(rgbImage, [192 82 60 52]);
subplot(2, 2, 2);
imshow(smallSubImage, []);
axis on;
title('Template Image to Search For');
% Search the red channel for a match.
correlationOutput = normxcorr2(smallSubImage(:,:,1), rgbImage(:,:,1));
subplot(2, 2, 3);
imshow(correlationOutput, []);
title('Correlation Output');
[maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
[ypeak, xpeak] = ind2sub(size(correlationOutput),maxIndex(1));
corr_offset = [(xpeak-size(smallSubImage,2)) (ypeak-size(smallSubImage,1))];
subplot(2, 2, 4);
imshow(rgbImage);
hold on;
rectangle('position',[corr_offset(1) corr_offset(2) 50 50],...
          'edgecolor','g','linewidth',2);
title('Template Image Found in Original Image');