classdef BodyAccelPDMixin
  properties
    robot
    use_mex
    numq
  end

  methods
    function obj = BodyAccelPDMixin(r, options)
      options = applyDefaults(options, ...
        struct('use_mex', 1));
      for f = fieldnames(options)'
        obj.(f{1}) = options.(f{1});
      end
      if (obj.use_mex && exist('bodyMotionControlmex','file')~=3)
          error('can''t find bodyMotionControlmex.  did you build it?');
      end

      obj.robot = r;
      obj.numq = obj.robot.getNumPositions();
    end

    function body_vdot = getBodyVdot(obj, t, x, body_data, params)

      if t < body_data.ts(1)
        t = body_data.ts(1);
      end
      if t > body_data.ts(2)
        t = body_data.ts(2);
      end
      [body_des, body_v_des, body_vdot_des] = evalCubicSplineSegment(t - body_data.ts(1), body_data.coefs);

      lcmgl = LCMGLClient(sprintf('link_%d_desired', body_data.body_id));
      lcmgl.sphere(body_des(1:3), 0.03, 20, 20);
      lcmgl.switchBuffers();

      if (obj.use_mex == 0 || obj.use_mex==2)
        q = x(1:obj.numq);
        qd = x(obj.numq+1:end);
        kinsol = doKinematics(obj.robot,q,false,true,qd);

        % TODO: this should be updated to use quaternions/spatial velocity
        [p,J] = forwardKin(obj.robot,kinsol,body_data.body_id,[0;0;0],1); 

        err = [body_des(1:3)-p(1:3);angleDiff(p(4:end),body_des(4:end))];

        body_vdot = params.Kp.*err + params.Kd.*(body_v_des-J*qd) + body_vdot_des; 
        if obj.use_mex == 2
          % check that matlab/mex agree
          body_vdot_mex = statelessBodyMotionControlmex(obj.robot.getMexModelPtr.ptr,x,body_data.body_id,body_des,body_v_des,body_vdot_des,params.Kp,params.Kd);  
          valuecheck(body_vdot_mex,body_vdot);
        end
      else
        body_vdot = statelessBodyMotionControlmex(obj.robot.getMexModelPtr.ptr,x,body_data.body_id,body_des,body_v_des,body_vdot_des,params.Kp,params.Kd);    
      end
    end
  end
end


