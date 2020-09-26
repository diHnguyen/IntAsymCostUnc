# IntAsymCostUnc
Network Interdiction with Asymmetric Cost Uncertainty - Di H. Nguyen and J. Cole Smith

Link to paper: [TBD]

Abstract: We study a network game between a leader and a follower having opposing objectives. The follower travels from source to sink using the least-cost path, while the leader tries to maximize the expected value of the evader’s cost. In this game, the leader only knows the cost distributions, not each arc’s exact value, and must use this information to make interdiction decisions in advance. An arc on the network carries an extra cost if attacked by the leader. The follower, after observing the attacks made on the network, has perfect information and traverses using a least-cost path. Selected decomposition-based algorithms and enhancements are discussed in our solution approach to this problem.  

In this repository you will find the following:

Test Instances:
- Preliminary Test Instances: Networks of 15 nodes with target density 0.2.
- Center Test Instances: Networks of 20 nodes with target density 0.2 + modified test instances based on Center Test Instances.

Algorithms and Enhancements:
- Algorithms: Base, Multicut, Multicut and Multipartition
- Enhancements: NoRep, Weighted, Tree-SA, Path-SA, Score 
(We do not apply Weighted enhancement to Multicut and Multipartition Algorithm.)
