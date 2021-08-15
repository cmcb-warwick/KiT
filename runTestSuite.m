% run a suite of tests to check basic functionality of KiT
% for tracking of Kinetochores in spinning disk or LLSM data
%
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%%%%%%%%%%%%%%%%%

%see here for details on unit tests and tags
%actually run via results = runtests({'runTestSuite'},'Tag','Diagnose');

classdef runTestSuite < matlab.unittest.TestCase
    methods (Test, TestTags = {'Data'})
        function loadExampleData (testCase)
            kitDataTest();
        end
    end
    methods (Test, TestTags = {'Detection'})
        function detection (testCase)
            kitExampleDetectionTest();
        end
        function refinement (testCase)
            [spots, movie, job] = kitExampleDetectionTest();
            kitExampleRefinementTest(spots,movie,job);
        end
    end
    methods (Test, TestTags = {'Plane'})
        function planeFit (testCase)
            kitExamplePlaneFitTest();
        end
    end
    methods (Test, TestTags = {'Tracking'})
        function tracking (testCase)
            kitExampleTrackTest(); 
        end
    end
    methods (Test, TestTags = {'Sisters'})
        function groupSisters (testCase)
            kitExampleGroupSistersTest();
        end
    end
    methods (Test, TestTags = {'Diagnose'})
        function diagnose (testCase)
            kitDiagnoseAndTest();
        end
    end
end
