options.floating = true;
N = 1;
p = TimeSteppingRigidBodyManipulator('Cylinder.urdf',.01,options);
for i=2:N,
  options.namesuffix = num2str(N);
  p = p.addRobotFromURDF('Cylinder.urdf',[],[],options);
end
p = p.addRobotFromURDF('ground_plane.urdf');

x0 = .1*randn(p.getNumDiscStates,1);
x0(3:6:end) = x0(3:6:end) + .5;

x0 = p.resolveConstraints(x0);
v = p.constructVisualizer();
v.drawWrapper(0,x0);

traj = simulate(p,[0 .5],x0);
v.playback(traj)
% function collisionDetectGradTest(visualize,n_debris)
%   if nargin < 1
%     visualize = false;
%   end
%   if nargin < 2
%     n_debris = 5;
%   end
%   n_debris = 5;
%   visualize = true;
%   options.floating = true;
%   S = warning('OFF','Drake:RigidBodyManipulator:UnsupportedContactPoints');
%   warning('OFF','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
%   warning('OFF','Drake:RigidBodyManipulator:UnsupportedJointLimits');
%   r = RigidBodyManipulator('Cylinder.urdf',options);
%   for i = 1:n_debris
%     r = r.addRobotFromURDF('FallingBrick.urdf',3*(2*rand(3,1)-1),2*pi*rand(3,1));
%     r = r.addRobotFromURDF('ball.urdf',3*(2*rand(3,1)-1),2*pi*rand(3,1));
%     r = r.addRobotFromURDF('Cylinder.urdf',3*(2*rand(3,1)-1),2*pi*rand(3,1));
%   end
%   r = r.addRobotFromURDF('ground_plane.urdf');
%   warning(S);
%   
%   if visualize
%     v = r.constructVisualizer();
%   end
%   
%   x0 = r.getInitialState;
%   x0 = r.resolveConstraints(x0);
% end