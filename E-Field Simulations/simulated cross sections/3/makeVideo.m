video = VideoWriter('Ez.avi'); %create the video object
video.FrameRate = 10;
open(video); %open the file for writing
for ii=1:239 %where N is the number of images
    sprintf('f%d.png',ii)
  I = imread(sprintf('f%d.png',ii)); %read the next image
  writeVideo(video,I); %write the image to file
end
close(video); %close the file