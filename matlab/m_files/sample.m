%% Sample Numerical Tour
% This tour defines and plot a discretized function.

%% Ploting a function

%%
% One can define the function
% \[f : x \in [-1,1] \rightarrow x^2 \in \mathbb{R}^+\]
% in matlab as a discretized vector.

f = linspace(-1,1,256).^2;

%EXO
%% As an exercise, display the function.
clf; h = plot(f); axis tight;
%EXO