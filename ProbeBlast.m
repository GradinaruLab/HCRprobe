function [num_hits, hits_pos, hits_names, bad_100hits] = ProbeBlast(seq, trgt_name)

seq_num = size(seq, 1);

if size(seq, 2) > 1
    seq_temp = [];
    
    for i=1:seq_num
        seqi = seq(i, :);
        seq_temp = [seq_temp; deblank(strcat(seqi(1), strcat(seqi(2), seqi(3))))];
    end
    seq = seq_temp;
end
seq_leng = strlength(seq(1));


%% making fasta file of the sequence list
% for i = 1:seq_num
%     seq_fa(i).Sequence = upper(seq(i));
%     seq_fa(i).Header = strcat(num2str(i));
% end
% fastawrite('seq_fasta.txt', seq_fa);


%% blasting the sequence list 
cmd = string('blastn -num_threads 10 -max_target_seqs 5000 -outfmt "6 qseqid stitle pident qcovs" ');

hits_pos = [];
hits_names = [];
bad_100hits = false(seq_num, 1);
num_hits = zeros(seq_num, 1);

mkdir blast_results

% seq_fa = struct([]);
disp('   (1) creating fasta files of each seq ...')
for i=1:seq_num
    seq_fa.Sequence = upper(seq(i));
    seq_fa.Header = strcat(num2str(i));
    
    fastawrite(sprintf('blast_results/seq_fasta_%d.txt', i), seq_fa);
end

disp('   (2) parallel blasting ...')
parfor i=1:seq_num
    if seq_leng < 30
        final_cmd = strcat(cmd, sprintf('-db mouse_transcript -task blastn-short -evalue 1000 -query blast_results/seq_fasta_%d.txt -out blast_results/blast_result_%d.txt', i, i))
        status = system(char(final_cmd));
    else
        final_cmd = strcat(cmd, sprintf('-db mouse_transcript -task blastn -query blast_results/seq_fasta_%d.txt -out blast_results/blast_result_%d.txt', i, i));
        status = system(char(final_cmd));
    end
end

disp('   (3) compiling results ...')
for i=1:seq_num
    %% reading blast result
    blast_result = importdata(sprintf('blast_results/blast_result_%d.txt', i));

    if ~isempty(blast_result)
        %% parsing results
        names = blast_result.textdata;
        values = blast_result.data;
        pident = values(:,1);
        qcovs = values(:,2);

        predicted = regexpi(names(:,2), 'predicted');
        predicted_inds = cellfun('isempty', predicted);
        hits_on = regexpi(names(:,2), trgt_name);
        hits_off_inds = cellfun('isempty', hits_on);
        off_inds = find(all([predicted_inds, hits_off_inds], 2));

        hits_name = names(off_inds, 2);
        pident = pident(off_inds);
        qcovs = qcovs(off_inds);
        hits_name_unique = unique(hits_name, 'stable');
        num_hits(i) = numel(hits_name_unique);

        hits_pos = [hits_pos; zeros(num_hits(i), 1) + i];
        hits_names = [hits_names; hits_name_unique];

        if ~isempty(find(pident==100 & qcovs==100))
            bad_100hits(i) = 1;
        end
    end
    
%     delete(sprintf('seq_fasta_%d.txt', i));
%     delete(sprintf('blast_result_%d.txt', i));
end
rmdir('blast_results', 's');
num_hits = num_hits';

% figure, plot(num_hits), hold on
% for i = 1:numel(bad_hits_idx)
%     if bad_hits_idx(i)
%         line([i,i], get('ylim'), 'linestyle', ':', 'color', 'r');
%     end
% end


