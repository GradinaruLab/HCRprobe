% function [prb_list, num_hits] = ProbeDesign(target, varargin)
tic

if isempty(gcp('nocreate'))
    parpool(20);
end

clearvars -except target

% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNA probe design for RNA ISH HCR v3 (split-probe)
%
% Default parameters
%
% probe length          :   20
% spacer length         :   2
% min space btwn probes :   3
% gc contents           :   45-60
%
% initiator number      :   1
% initiator type        :   I2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% input
% target = string('ATGGAGCGGAGACGCATCACCTCTGCGCGCCGCTCCTATGCCTCCGAGACGGTGGTCAGGGGCCTCGGTCCTAGTCGACAACTGGGTACCATGCCACGCTTCTCCTTGTCTCGAATGACTCCTCCACTCCCTGCCAGGGTGGACTTCTCCCTGGCCGGGGCGCTCAATGCTGGCTTCAAGGAGACACGGGCGAGCGAGCGTGCAGAGATGATGGAGCTCAATGACCGCTTTGCTAGCTACATCGAGAAGGTCCGCTTCCTGGAACAGCAAAACAAGGCGCTGGCAGCTGAACTGAACCAGCTTCGAGCCAAGGAGCCCACCAAACTGGCTGATGTCTACCAGGCGGAGCTTCGGGAGCTGCGGCTGCGGCTGGACCAGCTTACGGCCAACAGTGCCCGGCTGGAGGTGGAGAGGGACAACTTTGCACAGGACCTCGGCACCCTGAGGCAGAAGCTCCAAGATGAAACCAACCTGAGGCTGGAGGCAGAGAACAACCTGGCTGCGTATAGACAGGAGGCAGATGAAGCCACCCTGGCTCGTGTGGATTTGGAGAGAAAGGTTGAATCGCTGGAGGAGGAGATCCAGTTCTTAAGGAAGATCTATGAGGAGGAAGTTCGAGAACTCCGGGAGCAGCTGGCCCAACAGCAGGTCCACGTGGAGATGGATGTGGCCAAGCCAGACCTCACAGCGGCCCTGAGAGAGATTCGCACTCAATACGAGGCAGTGGCCACCAGTAACATGCAAGAGACAGAGGAGTGGTATCGGTCTAAGTTTGCAGACCTCACAGACGCTGCGTCCCGCAACGCAGAGCTGCTCCGCCAAGCCAAGCACGAAGCTAACGACTATCGCCGCCAACTGCAGGCCTTGACCTGCGATCTGGAGTCCCTGCGCGGCACGAACGAGTCCCTAGAGCGGCAAATGCGCGAACAGGAAGAGCGCCATGCGCGGGAGTCGGCCAGTTACCAGGAGGCACTTGCTCGGCTGGAGGAGGAGGGCCAAAGCCTCAAGGAGGAGATGGCCCGCCACCTGCAGGAGTACCAGGATCTACTCAACGTTAAGCTAGCCCTGGACATCGAGATCGCCACCTACAGGAAATTGCTGGAGGGCGAAGAAAACCGCATCACCATTCCTGTACAGACTTTCTCCAACCTCCAGATCCGAGAAACCAGCCTGGACACCAAATCCGTGTCAGAAGGCCACCTCAAGAGGAACATCGTGGTAAAGACTGTGGAGATGCGGGATGGTGAGGTCATTAAGGACTCGAAGCAGGAGCACAAGGACGTGGTGATGTGA');
trgt_name = 'calb2';
hairpin_idx = 4;


% initiator seqs
initseqs = {
    'gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg', 'gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg'; ...
    'CCTCgTAAATCCTCATCAATCATCCAgTAAACCgCC', 'AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC'; ...
    'gTCCCTgCCTCTATATCTCCACTCAACTTTAACCCg', 'AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC'; ...
    'CCTCAACCTACCTCCAACTCTCACCATATTCgCTTC', 'CACATTTACAgACCTCAACCTACCTCCAACTCTCAC'; ...
    'CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT', 'CACTTCATATCACTCACTCCCAATCTCTATCTACCC'};

spacers = {     % from Harry's paper (Development, 2016)
    'atat', ...
    'aaaa', ...
    'taaa', ...
    'aaaa', ...
    'attt'};

%% default params
prb_length = 20;
prb_length = prb_length*2;
spc_btwn = 3;
gc_range = [45,60];
initseq = char(initseqs(hairpin_idx, 2));
% spacer = spacers(hairpin_idx);
spacer = {'ta', 'at'} ;

% cdsonly = true;

%%%% sequence linearize
trgt = '';
for i=1:size(target, 1)
    trgt = strcat(trgt, target(i));
end


% %% get 
% gen_info = getgenbank(accessnum);
% seq = gen_info.Sequence;
% def = strsplit(string(gen_info.Definition), {'(', ')'});
% name = def(2);
% 
% if cdsonly
%     cds_info = gen_info.CDS;
%     cds_pos = cds_info.indices;
%     trgt = seq(cds_pos(1):cds_pos(2));
% else
%     trgt = seq;
% end


%% finding all potential probes
disp('- finding all potential probes ...');
comp_trgt = seqcomplement(char(trgt));

trgt_length = strlength(trgt);
index = repmat((0:trgt_length-prb_length)', 1, prb_length) + repmat(1:prb_length, trgt_length-prb_length+1, 1);
prb_list = seqreverse(comp_trgt(index));
prb_list_A = prb_list(:,1:prb_length/2);
prb_list_B = prb_list(:,prb_length/2+1:end);
prb_num = size(prb_list, 1);

%% (0-1) finding probes that have proper gc contents
disp(' 0-1. GC contents ...')
for i = trgt_length-prb_length+1:-1:1
    prb_props_A(i) = oligoprop(prb_list_A(i,:));
    prb_props_B(i) = oligoprop(prb_list_B(i,:));
end

gc_A = [prb_props_A.GC]';
gc_B = [prb_props_B.GC]';
bad_gc_A = gc_A < gc_range(1) | gc_A > gc_range(2);
bad_gc_B = gc_B < gc_range(1) | gc_B > gc_range(2);
bad_gc = bad_gc_A | bad_gc_B;


%% (0-2) finding probes with nucleotide repeats
disp(' 0-2. nucleotide repeats >= 4 ...')
repeats_A = regexpi(cellstr(prb_list_A), 'a{4,}|c{4,}|g{4,}|t{4,}', 'once');
repeats_B = regexpi(cellstr(prb_list_B), 'a{4,}|c{4,}|g{4,}|t{4,}', 'once');
bad_repeats_A = ~cellfun('isempty', repeats_A);
bad_repeats_B = ~cellfun('isempty', repeats_B);
bad_repeats = bad_repeats_A | bad_repeats_B;


%% (1) blasting target sequences
disp(' 1. blasting target sequences ...')
[num_hits, hits_pos, hits_names, bad_100hits] = ProbeBlast(string(prb_list), trgt_name);
figure, 
subplot(221), plot(num_hits, 'x-'); 
title('result of target sequence blast'), xlabel('probe no.'), ylabel('number of hits');


%% (2) finding free energy of secondary structure above the threshold
disp(' 2. dG of secondary structure ...')
dG_second_A = zeros(prb_num, 1);
dG_second_B = zeros(prb_num, 1);
parfor i=1:prb_num
    [est, dG_second_A(i)] = rnafold(prb_list_A(i,:));
    [est, dG_second_B(i)] = rnafold(prb_list_B(i,:));
end


%% attaching spacer and initiator
spcmat_A = repmat(spacer(2), prb_num, 1);
spcmat_B = repmat(spacer(1), prb_num, 1);
initmat_A = repmat(initseq(strlength(initseq)/2+1:end), prb_num, 1);
initmat_B = repmat(initseq(1:strlength(initseq)/2), prb_num, 1);
% prb_list_full_A = strcat(strcat(initmat_A, spcmat), prb_list_A);
% prb_list_full_B = strcat(prb_list_B, strcat(spcmat, initmat_B));
prb_list_full_A = strcat(prb_list_A, strcat(spcmat_A, initmat_A));
prb_list_full_B = strcat(strcat(initmat_B, spcmat_B), prb_list_B);


%% (3) blasting full sequences
disp(' 3. blasting full seq ...')
num_hits_full_A = ProbeBlast(string(prb_list_full_A), trgt_name);
num_hits_full_B = ProbeBlast(string(prb_list_full_B), trgt_name);
subplot(222), plot(num_hits_full_A, 'x-');
title('result of full sequence blast'), xlabel('probe no.'), ylabel('number of hits');


%% (4) checking cross-dimerization between probes
disp(' 4. cross-dimerization ...')
% [scr, prb_scr] = ProbeComp(prb_list);
[scr, prb_scr] = ProbeComp([char(prb_list_full_A); char(prb_list_full_B)]);
subplot(223), imagesc(scr), title('cross dimerization score'), colorbar, axis square


%% final decision
disp('- final decision!')
% thresholds
trgt_hits_thresh = 20;      % default:  20
dG_thresh = -9;             %           -9
offoverlap_thresh = 10;      %           5
full_hits_thresh = 20;      %           20
scr_thresh = (min(prb_scr)+max(prb_scr))/2;          %           10
% scr_thresh =10;

% from (1) blasting of target sequences
bad_hits = num_hits > trgt_hits_thresh;
bad_hits = bad_hits';

% from (2) min free energy of secondary structure
bad_second_A = dG_second_A < dG_thresh;
bad_second_B = dG_second_B < dG_thresh;
bad_second = bad_second_A | bad_second_B;

% from (3) blasting of full sequences
bad_hits_full_A = num_hits_full_A > full_hits_thresh;
bad_hits_full_B = num_hits_full_B > full_hits_thresh;
bad_hits_full_A = bad_hits_full_A';
bad_hits_full_B = bad_hits_full_B';
bad_hits_full = bad_hits_full_A | bad_hits_full_B;

% from (4) cross-dimerization check
bad_scr_A = prb_scr(1:size(prb_scr,1)/2) > scr_thresh;
bad_scr_B = prb_scr(size(prb_scr,1)/2+1:end) > scr_thresh;
bad_scr = bad_scr_A | bad_scr_B;
% bad_scr = false(prb_num, 1);
% 
% bad_repeats = false(prb_num, 1);
% bad_second = false(prb_num, 1);
bad_inds = [bad_gc, bad_repeats, bad_100hits, bad_hits, bad_second, bad_hits_full, bad_scr];
num_bads = sum(bad_inds, 1)
prb_pos = find(all(~bad_inds, 2));
if isempty(prb_pos)
    disp('No probe detected!')
    return;
end

prb_hits_names = [];
prb_hits_pos = [];
for i=1:numel(prb_pos)
    hits_names_i = hits_names(find(hits_pos == prb_pos(i)));
    prb_hits_names = [prb_hits_names; hits_names_i];
    prb_hits_pos = [prb_hits_pos; zeros(numel(hits_names_i), 1)+prb_pos(i)];
end

max_num_offtarget_rec = [];
while 1
    [hits_names_uniq, ia, ic] = unique(prb_hits_names, 'stable');
    offoverlap_num = histcounts(ic, numel(ia));
    [max_num_offtarget, max_num_offtarget_ind] = max(offoverlap_num);
    max_num_offtarget_rec = [max_num_offtarget_rec; max_num_offtarget];
    subplot(224), plot(max_num_offtarget_rec, 'x-'),
    xlabel('iteration'), ylabel('max num of probes targeting off'),
    title('probe dropping')
    
    if max_num_offtarget <= offoverlap_thresh
        break;
    end
    
    prb_pos_drop = prb_hits_pos(find(ic == max_num_offtarget_ind(1)));
    drop_ind = [];
    for i=1:numel(prb_pos_drop)
        drop_ind = [drop_ind; find(prb_hits_pos == prb_pos_drop(i))];
    end
    prb_hits_names(drop_ind) = [];
    prb_hits_pos(drop_ind) = [];
end

prb_pos_sorted = unique(prb_hits_pos);
prb_hits_names_uniq = unique(prb_hits_names);
prb_final_pos = prb_pos_sorted(1);
for i=2:numel(prb_pos_sorted)
    distance = abs(prb_final_pos - prb_pos_sorted(i));
    if isempty(find(distance < (prb_length + spc_btwn), 1))
        prb_final_pos = [prb_final_pos; prb_pos_sorted(i)];
    end
end

prb_hits_names_selected = [];
for i=1:numel(prb_final_pos)
    hits_names_i_selected = hits_names(find(hits_pos == prb_final_pos(i)));
    prb_hits_names_selected = [prb_hits_names_selected; hits_names_i_selected];
end
[hits_names_selected_uniq, ia, ic] = unique(prb_hits_names_selected, 'stable');
offoverlap_num_selected = histcounts(ic, numel(ia));
max_num_offtarget_final = max(offoverlap_num_selected);

max_num_offtarget_rec
max_num_offtarget_final
prb_final_num = numel(prb_final_pos)
prb_final_seq = [repmat(initseq(1:strlength(initseq)/2), prb_final_num, 1), ...
    repmat(spacer(1), prb_final_num, 1), ...
    string(prb_list_B(prb_final_pos, :)), ...
    string(prb_list_A(prb_final_pos, :)), ...
    repmat(spacer(2), prb_final_num, 1), ...
    repmat(initseq(strlength(initseq)/2+1:end), prb_final_num, 1)];
% prb_final_seq = [string(prb_list(prb_final_pos, :)), repmat(spacer, prb_final_num, 1), repmat(initseq, prb_final_num, 1)];


toc
