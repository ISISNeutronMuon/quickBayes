Welcome To quickBayes
=====================

The quickBayes package is an open source library for calculating Bayesian quantities in a short period of time, by making some assumptions.
The package is cross platform, supporting Windows, Mac OS and Linux.
This package has been developed by Anthony Lim from STFC’s ISIS Neutron and Muon facility.

The quickBayes package was originally designed as a successor to the quasielasticbayes package.
However, quickBayes has abstracted the key ideas to make a more generic package.

In science hypotheses are tested against data, to determine the underlying behavior of the system.
These hypotheses can be in the form of a mathematical expression that originates from first principles (i.e. it has been derived) or is an approximation to other more complex mechanisms (e.g. semi-empirical methods).
The quickBayes package is designed to test these hypotheses against the user’s data to identify the most likely.
For example, in Quasi Elastic Neutron Scattering (QENS) the data can be represented by a summation of Lorentzian peaks (convoluted with a resolution function, see  \textbf{link to doc} for more detail), but the number of peaks is not obvious as shown by \textbf{insert figure}.
The quickBayes packages is designed to make these decisions easier by calculating the probability of one, two or three peaks given the data.

This documentation is split into four main parts:
\begin{enumerate}
\item An introduction to Bayesian Statistics \textbf{link}
\item The key principles behind quickBayes \textbf{link}
\item Some real world examples \textbf{link} from QENS and muon spectroscopy
\item Developer documentation \textbf{link}
\end{enumerate}
\end{abstract}
