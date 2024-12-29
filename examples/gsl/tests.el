(require 'cl-macs)
(add-to-list 'load-path (expand-file-name "."))


(defun test-linalg () 
  (require 'gsl-linalg)

  (let ((ans (gsl-linalg-LU-solve [[1.0 2.0]
				   [3. 4.]]
				  [23.0 53.0])))
    ;; (cl-assert
    ;;  (equal ans [7.0 8.0]))
    (print (format "You should get [7.0 8.0]"))
    (print (format "You actually got %s" ans))))

(defun test-dgemm ()
  (require 'gsl-linalg)

  (let ((result (gsl-blas-dgemm [[0.18 0.60 0.57 0.96]
				 [0.41 0.24 0.99 0.58]
				 [0.14 0.30 0.97 0.66]
				 [0.51 0.13 0.19 0.85]]
				[[-4.052050229573973]
				 [-12.605611395906903]
				 [1.6609116267088417]
				 [8.693766928795227]])))

    (print result)
    (print "You should get this:\n[[1.0000000000000004] [1.9999999999999991] [2.999999999999999] [4.0]]")
    (print "No error counts as a passed test here.")))



