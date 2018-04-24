function [scr, prb_scr] = ProbeComp(seq)

% seq = [seq; 'aaaaaaaaaaaaaaaaaaaa', 'aaaa', 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'];
% seq = seq(25:42,:);
scr_mat = [-1, -1, -1, 1; -1, -1, 1, -1; -1, 1, -1, -1; 1, -1, -1, -1];
scr = zeros(size(seq,1));
% 
% max_scr = 0;
% for i= 1:size(seq, 1)
%     for j = 1:size(seq, 1)
%         if i<=j
%             seqi = seq(i,:);
%             seqi = deblank(strcat(seqi(1), strcat(seqi(2), seqi(3))));
%             seqj = seq(j,:);
%             seqj = deblank(strcat(seqj(1), strcat(seqj(2), seqj(3))));
% 
%     %         scr(i,j) = swalign(char(seqi), char(seqj), ...
%     %             'scoringmatrix', scr_mat, 'gapopen', 5, 'alpha', 'nt');
% 
%             scr(i,j) = swalign(char(seqi), char(seqj), 'alpha', 'nt');
%             if max_scr < scr(i,j)
%                 max_scr = scr(i,j);
%             end
%         else
%             scr(i,j) = max_scr;
%         end
%             
%     end
% end
% 
% figure, imagesc(scr), title('cross dimerization score'), colorbar


for i=1:size(seq, 1)
    seqi = seq(i,:);
    parfor j=1:size(seq, 1)
        seqj = seq(j,:);
        scr(i,j) = swalign(char(seqi), seqrcomplement(char(seqj)), ...
            'scoringmatrix', scr_mat, 'gapopen', 5, 'alpha', 'nt');
%         scr(i,j) = swalign(char(seqi), seqrcomplement(char(seqj)), ...
%             'gapopen', 5, 'alpha', 'nt');
    end
end
prb_scr = max(scr, [], 2);
            