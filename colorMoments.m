function colorMoments = colorMoments(image)
% input: image to be analyzed and extract 2 first moments from each R,G,B
% output: 1x6 vector containing the 2 first color momenst from each R,G,B
% channel

% extract color channels
R = double(image(:, :, 1));
G = double(image(:, :, 2));
B = double(image(:, :, 3));

% compute 2 first color moments from each channel
meanR = mean( R(:) );
stdR  = std( R(:) ); 
skewR = skewness( R(:) );
meanG = mean( G(:) );
stdG  = std( G(:) );
skewG = skewness( G(:) );
meanB = mean( B(:) );
stdB  = std( B(:) );
skewB = skewness( B(:) );

% construct output vector
colorMoments = zeros(1, 9);
colorMoments(1, :) = [meanR stdR skewR meanG stdG skewG meanB stdB skewB];

% clear workspace
clear('R', 'G', 'B', 'meanR', 'stdR', 'meanG', 'stdG', 'meanB', 'stdB');

end