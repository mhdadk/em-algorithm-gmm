function frames = save_GIF(X,num_iter,mu_save,Sigma_save,comp_PMF_save)

% create figure handle

h = figure;

% initialize structure to store animation frames

frames(num_iter) = struct('cdata',[],'colormap',[]);

% plot vectors sampled from GMM

plot(X(:,1),X(:,2),'+')
set(gcf,'color','w')
lim = axis;
axis([lim(1)-1,lim(2)+1,lim(3)-1,lim(4)+1])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')

% ensure that getframe() returns a consistent size

axis manual

% name of .gif file to be saved

filename = 'GMM.gif';

% suppress optimization suggestion warning from MATLAB

warning('off','MATLAB:fplot:NotVectorized')

for i = 1:num_iter
    
    % create a gmdistribution object and get its PDF
    
    gm = gmdistribution(mu_save(:,:,i),...
                        Sigma_save(:,:,:,i),...
                        comp_PMF_save(:,:,i));
                    
    gm_PDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
    
    % plot the contours for the PDF of the GMM on the same figure
    
    hold on
    fc = fcontour(h.CurrentAxes,gm_PDF);
    hold off
    
    % capture the plot as an image
    
    frames(i) = getframe(h);
    
    % delete the current contour plot to make way for the next one in the
    % next iteration
    
    delete(fc);
    
    % convert frame to image
    
    im = frame2im(frames(i)); 
    
    % convert image to indexed image to save as .gif file
    
    [imind,cm] = rgb2ind(im,256);
    
    % write to the .gif File
    
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end

% close the figure to make way for the animated plot

close(h)

end