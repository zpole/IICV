function experiments()
% EXPERIMENTS   Run image classification experiments
%    The experimens download a number of benchmark datasets in the
%    'data/' subfolder. Make sure that there are several GBs of
%    space available.
%
%    By default, experiments run with a lite option turned on. This
%    quickly runs all of them on tiny subsets of the actual data.
%    This is used only for testing; to run the actual experiments,
%    set the lite variable to false.
%
%    Running all the experiments is a slow process. Using parallel
%    MATLAB and several cores/machiens is suggested.

% Author: Andrea Vedaldi

% Copyright (C) 2013 Andrea Vedaldi
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

lite = true ;
clear ex ;

ex(1).prefix = 'fv-aug' ;
ex(1).trainOpts = {'C', 10} ;
ex(1).datasets = {'fmd', 'scene67'} ;
ex(1).seed = 1 ;
ex(1).opts = {...
  'type', 'fv', ...
  'numWords', 256, ...
  'layouts', {'1x1'}, ...
  'geometricExtension', 'xy', ...
  'numPcaDimensions', 80, ...
  'extractorFn', @(x) getDenseSIFT(x, ...
                                   'step', 4, ...
                                   'scales', 2.^(1:-.5:-3))};

ex(2) = ex(1) ;
ex(2).datasets = {'caltech101'}  ;
ex(2).opts{end} =  @(x) getDenseSIFT(x, ...
                                     'step', 4, ...
                                     'scales', 2.^(0:-.5:-3)) ;

ex(3) = ex(1) ;
ex(3).datasets = {'voc07'}  ;
ex(3).C = 1 ;

ex(4) = ex(1) ;
ex(4).prefix = 'vlad-aug' ;
ex(4).opts = {...
  'type', 'vlad', ...
  'numWords', 256, ...
  'layouts', {'1x1'}, ...
  'geometricExtension', 'xy', ...
  'numPcaDimensions', 100, ...
  'whitening', true, ...
  'whiteningRegul', 0.01, ...
  'renormalize', true, ...
  'extractorFn', @(x) getDenseSIFT(x, ...
                                   'step', 4, ...
                                   'scales', 2.^(1:-.5:-3))};
ex(5) = ex(4) ;
ex(5).datasets = {'caltech101'} ;
ex(5).opts{end} = ex(2).opts{end} ;

ex(6) = ex(4) ;
ex(6).datasets = {'voc07'} ;
ex(6).C = 1 ;

ex(7) = ex(1) ;
ex(7).prefix = 'bovw-aug' ;
ex(7).opts = {...
  'type', 'bovw', ...
  'numWords', 4096, ...
  'layouts', {'1x1'}, ...
  'geometricExtension', 'xy', ...
  'numPcaDimensions', 100, ...
  'whitening', true, ...
  'whiteningRegul', 0.01, ...
  'renormalize', true, ...
  'extractorFn', @(x) getDenseSIFT(x, ...
                                   'step', 4, ...
                                   'scales', 2.^(1:-.5:-3))};

ex(8) = ex(7) ;
ex(8).datasets = {'caltech101'} ;
ex(8).opts{end} = ex(2).opts{end} ;

ex(9) = ex(7) ;
ex(9).datasets = {'voc07'} ;
ex(9).C = 1 ;

ex(10).prefix = 'fv' ;
ex(10).trainOpts = {'C', 10} ;
ex(10).datasets = {'fmd', 'scene67'} ;
ex(10).seed = 1 ;
ex(10).opts = {...
  'type', 'fv', ...
  'numWords', 256, ...
  'layouts', {'1x1'}, ...
  'geometricExtension', 'none', ...
  'numPcaDimensions', 80, ...
  'extractorFn', @(x) getDenseSIFT(x, ...
                                   'step', 4, ...
                                   'scales', 2.^(1:-.5:-3))};

ex(11) = ex(10) ;
ex(11).datasets = {'caltech101'}  ;
ex(11).opts{end} =  @(x) getDenseSIFT(x, ...
                                     'step', 4, ...
                                     'scales', 2.^(0:-.5:-3)) ;

ex(12) = ex(10) ;
ex(12).datasets = {'voc07'}  ;
ex(12).C = 1 ;


ex(13).prefix = 'fv-sp' ;
ex(13).trainOpts = {'C', 10} ;
ex(13).datasets = {'fmd', 'scene67'} ;
ex(13).seed = 1 ;
ex(13).opts = {...
  'type', 'fv', ...
  'numWords', 256, ...
  'layouts', {'1x1', '3x1'}, ...
  'geometricExtension', 'none', ...
  'numPcaDimensions', 80, ...
  'extractorFn', @(x) getDenseSIFT(x, ...
                                   'step', 4, ...
                                   'scales', 2.^(1:-.5:-3))};

ex(14) = ex(13) ;
ex(14).datasets = {'caltech101'}  ;
ex(14).opts{6} = {'1x1', '2x2'} ;
ex(14).opts{end} =  @(x) getDenseSIFT(x, ...
                                     'step', 4, ...
                                     'scales', 2.^(0:-.5:-3)) ;

ex(15) = ex(13) ;
ex(15).datasets = {'voc07'}  ;
ex(15).C = 1 ;


if lite, tag = 'lite' ;
else tag = 'ex' ; end

i=6;
%j=1;
%for i=1:numel(ex)
for j=1:numel(ex(i).datasets)
    dataset = ex(i).datasets{j} ;
    if ~isfield(ex(i), 'trainOpts') || ~iscell(ex(i).trainOpts)
      ex(i).trainOpts = {} ;
    end
    traintest(...
      'prefix', [tag '-' dataset '-' ex(i).prefix], ...
      'seed', ex(i).seed, ...
      'dataset', char(dataset), ...
      'datasetDir', fullfile('data', dataset), ...
      'lite', lite, ...
      ex(i).trainOpts{:}, ...
      'encoderParams', ex(i).opts) ;
%  end
end

% print HTML table
%pf('<table>\n') ;
%ph('method', 'VOC07', 'Caltech 101', 'Scene 67', 'FMD') ;
%pr('FV', ...
%   ge([tag '-voc07-fv'],'ap11'), ...
%   ge([tag '-caltech101-fv']), ...
%   ge([tag '-scene67-fv']), ...
%   ge([tag '-fmd-fv'])) ;
%pr('FV + aug.', ...
%   ge([tag '-voc07-fv-aug'],'ap11'), ...
%   ge([tag '-caltech101-fv-aug']), ...
%   ge([tag '-scene67-fv-aug']), ...
%   ge([tag '-fmd-fv-aug'])) ;
%pr('FV + s.p.', ...
%   ge([tag '-voc07-fv-sp'],'ap11'), ...
%   ge([tag '-caltech101-fv-sp']), ...
%   ge([tag '-scene67-fv-sp']), ...
%   ge([tag '-fmd-fv-sp'])) ;
%pr('VLAD', ...
%   ge([tag '-voc07-vlad'],'ap11')), ...
%   ge([tag '-caltech101-vlad']), ...
%   ge([tag '-scene67-vlad']), ...
%   ge([tag '-fmd-vlad'])) ;
%pr('VLAD + aug.', ...
%   ge([tag '-voc07-vlad-aug'],'ap11'))%, ...
%   ge([tag '-caltech101-vlad-aug']), ...
%   ge([tag '-scene67-vlad-aug']), ...
%   ge([tag '-fmd-vlad-aug'])) ;
%pr('VLAD+sp', ...
%   ge([tag '-voc07-vlad-sp'],'ap11'), ...
%   ge([tag '-caltech101-vlad-sp']), ...
%   ge([tag '-scene67-vlad-sp']), ...
%   ge([tag '-fmd-vlad-sp'])) ;
%pr('BOVW', ...
%   ge([tag '-voc07-bovw'],'ap11'), ...
%   ge([tag '-caltech101-bovw']), ...
%   ge([tag '-scene67-bovw']), ...
%   ge([tag '-fmd-bovw'])) ;
%pr('BOVW + aug.', ...
%  ge([tag '-voc07-bovw-aug'],'ap11'), ...
%  ge([tag '-caltech101-bovw-aug']), ...
%  ge([tag '-scene67-bovw-aug']), ...
%  ge([tag '-fmd-bovw-aug'])) ;
%pr('BOVW+sp', ...
%   ge([tag '-voc07-bovw-sp'],'ap11'), ...
%   ge([tag '-caltech101-bovw-sp']), ...
%   ge([tag '-scene67-bovw-sp']), ...
%   ge([tag '-fmd-bovw-sp'])) ;
%pf('</table>\n');

function pf(str)
fprintf(str) ;

function str = ge(name, format)
if nargin == 1, format = 'acc'; end
data = load(fullfile('data', name, 'result.mat')) ;
switch format
  case 'acc'
    str = sprintf('%.2f%% <span style="font-size:8px;">Acc</span>', mean(diag(data.confusion)) * 100) ;
  case 'ap11'
    str = sprintf('%.2f%% <span style="font-size:8px;">mAP</span>', mean(data.ap11) * 100) ;
end

function pr(varargin)
fprintf('<tr>') ;
for i=1:numel(varargin), fprintf('<td>%s</td>',varargin{i}) ; end
fprintf('</tr>\n') ;

function ph(varargin)
fprintf('<tr>') ;
for i=1:numel(varargin), fprintf('<th>%s</th>',varargin{i}) ; end
fprintf('</tr>\n') ;
