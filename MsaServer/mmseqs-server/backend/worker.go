package main

import (
	"archive/tar"
	"bufio"
	"compress/gzip"
	"encoding/json"
	"errors"
	"fmt"
	"log"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"strconv"
	"strings"
	"sync/atomic"
	"syscall"
	"time"
)

type JobExecutionError struct {
	internal error
}

func (e *JobExecutionError) Error() string {
	return "Execution Error: " + e.internal.Error()
}

type JobTimeoutError struct {
}

func (e *JobTimeoutError) Error() string {
	return "Timeout"
}

type JobInvalidError struct {
}

func (e *JobInvalidError) Error() string {
	return "Invalid"
}

func execCommand(verbose bool, parameters ...string) (*exec.Cmd, chan error, error) {
	cmd := exec.Command(
		parameters[0],
		parameters[1:]...,
	)

	SetSysProcAttr(cmd)

	// Make sure MMseqs2's progress bar doesn't break
	cmd.Env = append(os.Environ(), "TTY=0")

	if verbose {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
	}

	done := make(chan error, 1)
	err := cmd.Start()
	if err != nil {
		return cmd, done, err
	}

	go func() {
		done <- cmd.Wait()
	}()

	return cmd, done, err
}

func RunJob(request JobRequest, config ConfigRoot) (err error) {
	switch job := request.Job.(type) {
	case SearchJob:
		resultBase := filepath.Join(config.Paths.Results, string(request.Id))
		for _, database := range job.Database {
			params, err := ReadParams(filepath.Join(config.Paths.Databases, database+".params"))
			if err != nil {
				return &JobExecutionError{err}
			}
			columns := "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln"
			if params.Taxonomy {
				columns += ",taxid,taxname"
			}
			parameters := []string{
				config.Paths.Mmseqs,
				"easy-search",
				filepath.Join(resultBase, "job.fasta"),
				filepath.Join(config.Paths.Databases, database),
				filepath.Join(resultBase, "alis_"+database),
				filepath.Join(resultBase, "tmp"),
				"--shuffle",
				"0",
				"--db-output",
				"--db-load-mode",
				"2",
				"--write-lookup",
				"1",
				"--format-output",
				columns,
			}
			parameters = append(parameters, strings.Fields(params.Search)...)

			if job.Mode == "summary" {
				parameters = append(parameters, "--greedy-best-hits")
			}

			if params.Taxonomy && job.TaxFilter != "" {
				parameters = append(parameters, "--taxon-list")
				parameters = append(parameters, job.TaxFilter)
			}

			cmd, done, err := execCommand(config.Verbose, parameters...)
			if err != nil {
				return &JobExecutionError{err}
			}

			select {
			case <-time.After(1 * time.Hour):
				if err := KillCommand(cmd); err != nil {
					log.Printf("Failed to kill: %s\n", err)
				}
				return &JobTimeoutError{}
			case err := <-done:
				if err != nil {
					return &JobExecutionError{err}
				}
			}
		}

		path := filepath.Join(filepath.Clean(config.Paths.Results), string(request.Id))
		file, err := os.Create(filepath.Join(path, "mmseqs_results_"+string(request.Id)+".tar.gz"))
		if err != nil {
			return &JobExecutionError{err}
		}
		err = ResultArchive(file, request.Id, path)
		if err != nil {
			file.Close()
			return &JobExecutionError{err}
		}
		err = file.Close()
		if err != nil {
			return &JobExecutionError{err}
		}

		if config.Verbose {
			log.Print("Process finished gracefully without error")
		}
		return nil
	case StructureSearchJob:
		resultBase := filepath.Join(config.Paths.Results, string(request.Id))
		for _, database := range job.Database {
			params, err := ReadParams(filepath.Join(config.Paths.Databases, database+".params"))
			if err != nil {
				return &JobExecutionError{err}
			}
			var mode2num = map[string]string{"3di": "0", "tmalign": "1", "3diaa": "2"}
			mode, found := mode2num[job.Mode]
			if !found {
				return &JobExecutionError{errors.New("Invalid mode selected")}
			}
			columns := "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln,tca,tseq"
			if params.Taxonomy {
				columns += ",taxid,taxname"
			}
			parameters := []string{
				config.Paths.FoldSeek,
				"easy-search",
				filepath.Join(resultBase, "job.pdb"),
				filepath.Join(config.Paths.Databases, database),
				filepath.Join(resultBase, "alis_"+database),
				filepath.Join(resultBase, "tmp"),
				// "--shuffle",
				// "0",
				"--alignment-type",
				mode,
				"--db-output",
				"--db-load-mode",
				"2",
				"--write-lookup",
				"1",
				"--format-output",
				columns,
			}
			parameters = append(parameters, strings.Fields(params.Search)...)

			if job.Mode == "summary" {
				parameters = append(parameters, "--greedy-best-hits")
			}

			if params.Taxonomy && job.TaxFilter != "" {
				parameters = append(parameters, "--taxon-list")
				parameters = append(parameters, job.TaxFilter)
			}

			cmd, done, err := execCommand(config.Verbose, parameters...)
			if err != nil {
				return &JobExecutionError{err}
			}

			select {
			case <-time.After(1 * time.Hour):
				if err := KillCommand(cmd); err != nil {
					log.Printf("Failed to kill: %s\n", err)
				}
				return &JobTimeoutError{}
			case err := <-done:
				if err != nil {
					return &JobExecutionError{err}
				}
			}
		}

		path := filepath.Join(filepath.Clean(config.Paths.Results), string(request.Id))
		file, err := os.Create(filepath.Join(path, "mmseqs_results_"+string(request.Id)+".tar.gz"))
		if err != nil {
			return &JobExecutionError{err}
		}
		err = ResultArchive(file, request.Id, path)
		if err != nil {
			file.Close()
			return &JobExecutionError{err}
		}
		err = file.Close()
		if err != nil {
			return &JobExecutionError{err}
		}

		if config.Verbose {
			log.Print("Process finished gracefully without error")
		}
		return nil
	case MsaJob:
		resultBase := filepath.Join(config.Paths.Results, string(request.Id))

		scriptPath := filepath.Join(resultBase, "msa.sh")
		script, err := os.Create(scriptPath)
		if err != nil {
			return &JobExecutionError{err}
		}

		script.WriteString(`#!/bin/bash -e
	MMSEQS="$1"
	QUERY="$2"
	BASE="$4"
	DB1="$5"
	DB2="$6"
	DB3="$7"
	USE_ENV="$8"
	USE_TEMPLATES="$9"
	FILTER="${10}"
	TAXONOMY="${11}"
	M8OUT="${12}"
	OUT="$13"
	
	mkdir -p "${BASE}"
	SEARCH_PARAM="--num-iterations 3 --db-load-mode 2 -a --k-score 'seq:96,prof:80' -e 0.1 --max-seqs 10000"
	EXPAND_PARAM="--expansion-mode 0 -e inf --expand-filter-clusters 0 --max-seq-id 0.95"
	export MMSEQS_CALL_DEPTH=1
	"${MMSEQS}" createdb "${QUERY}" "${BASE}/qdb" --shuffle 0
	python3 mmseqs-server/backend/aln_or_a3mtax.py "${BASE}/job.fasta"

		if [! -f "${BASE}/ALN_FOUND"]; then
			echo CALCULATING ALN
			"${MMSEQS}" search "${BASE}/qdb" "${DB1}" "${BASE}/res" "${BASE}/tmp" $SEARCH_PARAM
			"${MMSEQS}" expandaln "${BASE}/qdb" "${DB1}.idx" "${BASE}/res" "${DB1}.idx" "${BASE}/res_exp" --db-load-mode 2 ${EXPAND_PARAM}
			"${MMSEQS}" align   "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp" "${BASE}/res_exp_realign" --db-load-mode 2 -e 0.001 --max-accept 1000000 -c 0.5 --cov-mode 1
			"${MMSEQS}" cpdb "${BASE}/qdb.lookup" "${BASE}/res_exp_realign.lookup"
			"${MMSEQS}" unpackdb "${BASE}/res_exp_realign" "${BASE}" --unpack-name-mode 1 --unpack-suffix .aln
			"${MMSEQS}" rmdb "${BASE}/qdb"
			"${MMSEQS}" rmdb "${BASE}/qdb_h"
			"${MMSEQS}" rmdb "${BASE}/res"
			"${MMSEQS}" rmdb "${BASE}/res_exp"
			"${MMSEQS}" rmdb "${BASE}/res_final"
			"${MMSEQS}" rmdb "${BASE}/res_exp_realign"
			rm -rf -- "${BASE}/tmp"
			cd "${BASE}"
			tar -czvf "mmseqs_results_${OUT}.tar.gz" *.aln msa.sh
		else 
			echo CALCULATING A3M FILES
			"${MMSEQS}" convertalis "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign" "${BASE}/convertalis_tax" --format-output target,evalue,taxid,taxname,taxlineage --db-load-mode 2
			"${MMSEQS}" convertalis "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign" "${BASE}/convertalis_seq" --format-output target,tseq --db-load-mode 2
			"${MMSEQS}" result2msa "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign" "${BASE}/uniref.a3m" --msa-format-mode 6 --db-load-mode 2
			"${MMSEQS}" mvdb "${BASE}/uniref.a3m" "${BASE}/final.a3m"
			"${MMSEQS}" cpdb "${BASE}/qdb.lookup" "${BASE}/final.a3m.lookup" 
			"${MMSEQS}" unpackdb "${BASE}/final.a3m" "${BASE}" --unpack-name-mode 1 --unpack-suffix .a3m
			"${MMSEQS}" rmdb "${BASE}/final.a3m"
			python3 mmseqs-server/backend/add_tax_to_msa.py "${BASE}/convertalis_tax" "${BASE}"
			python3 mmseqs-server/backend/convertalis_seq_to_tsv.py "${BASE}/convertalis_seq" --a3m_dir "${BASE}"
			"${MMSEQS}" rmdb "${BASE}/qdb"
			"${MMSEQS}" rmdb "${BASE}/qdb_h"
			"${MMSEQS}" rmdb "${BASE}/res_exp_realign"
			"${MMSEQS}" rmdb "${BASE}/convertalis_tax"
			"${MMSEQS}" rmdb "${BASE}/convertalis_seq"
			rm -rf -- "${BASE}/tmp"
			cd "${BASE}"
			tar -czvf "mmseqs_results_${OUT}.tar.gz" *.a3m *.a3m.tax *.tsv msa.sh
		fi
		`)		


		// if [ ! -f "${BASE}/ALN_FOUND" ]; then
		// echo CALCULATING ALN
		// "${MMSEQS}" search "${BASE}/qdb" "${DB1}" "${BASE}/res" "${BASE}/tmp" $SEARCH_PARAM
		// "${MMSEQS}" expandaln "${BASE}/qdb" "${DB1}.idx" "${BASE}/res" "${DB1}.idx" "${BASE}/res_exp" --db-load-mode 2 ${EXPAND_PARAM}
		// "${MMSEQS}" align   "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp" "${BASE}/res_exp_realign" --db-load-mode 2 -e 0.001 --max-accept 1000000 -c 0.5 --cov-mode 1
		// "${MMSEQS}" cpdb "${BASE}/qdb.lookup" "${BASE}/res_exp_realign.lookup"
		// "${MMSEQS}" unpackdb "${BASE}/res_exp_realign" "${BASE}" --unpack-name-mode 1 --unpack-suffix .aln
		// "${MMSEQS}" rmdb "${BASE}/qdb"
		// "${MMSEQS}" rmdb "${BASE}/qdb_h"
		// "${MMSEQS}" rmdb "${BASE}/res"
		// "${MMSEQS}" rmdb "${BASE}/res_exp"
		// "${MMSEQS}" rmdb "${BASE}/res_final"
		// "${MMSEQS}" rmdb "${BASE}/res_exp_realign"
		// rm -rf -- "${BASE}/tmp"
		// cd "${BASE}"
		// tar -czvf "mmseqs_results_${OUT}.tar.gz" *.aln msa.sh

		err = script.Close()
		if err != nil {
			return &JobExecutionError{err}
		}

		modes := strings.Split(job.Mode, "-")
		useEnv := isIn("env", modes) != -1
		useTemplates := isIn("notemplates", modes) == -1
		useFilter := isIn("nofilter", modes) == -1
		taxonomy := isIn("taxonomy", modes) == 1
		m8out := isIn("m8output", modes) == 1
		var b2i = map[bool]int{false: 0, true: 1}

		parameters := []string{
			"/bin/sh",
			scriptPath,
			config.Paths.Mmseqs,
			filepath.Join(resultBase, "job.fasta"),
			"",
			resultBase,
			config.Paths.ColabFold.Uniref,
			config.Paths.ColabFold.Pdb,
			config.Paths.ColabFold.Environmental,
			strconv.Itoa(b2i[useEnv]),
			strconv.Itoa(b2i[useTemplates]),
			strconv.Itoa(b2i[useFilter]),
			strconv.Itoa(b2i[taxonomy]),
			strconv.Itoa(b2i[m8out]),
			string(request.Id),
			resultBase,
		}

		cmd, done, err := execCommand(config.Verbose, parameters...)
		if err != nil {
			return &JobExecutionError{err}
		}

		select {
		case <-time.After(1 * time.Hour):
			if err := KillCommand(cmd); err != nil {
				log.Printf("Failed to kill: %s\n", err)
			}
			return &JobTimeoutError{}
		case err := <-done:
			if err != nil {
				return &JobExecutionError{err}
			}

			path := filepath.Join(filepath.Clean(config.Paths.Results), string(request.Id))
			file, err := os.Create(filepath.Join(path, "ignore_mmseqs_results_"+string(request.Id)+".tar.gz"))
			if err != nil {
				return &JobExecutionError{err}
			}

			err = func() (err error) {
				gw := gzip.NewWriter(file)
				defer func() {
					cerr := gw.Close()
					if err == nil {
						err = cerr
					}
				}()
				tw := tar.NewWriter(gw)
				defer func() {
					cerr := tw.Close()
					if err == nil {
						err = cerr
					}
				}()

				/*
				if config.App == AppPredictProtein {
					if err := addFile(tw, filepath.Join(resultBase, "uniref.sto")); err != nil {
						return err
					}

					if err := addFile(tw, filepath.Join(resultBase, "uniref.m8")); err != nil {
						return err
					}

					if err := addFile(tw, filepath.Join(resultBase, "pdb70.sto")); err != nil {
						return err
					}

					if err := addFile(tw, filepath.Join(resultBase, "pdb70.m8")); err != nil {
						return err
					}
				} else {
					suffix := ".a3m"
					if m8out {
						suffix = ".m8"
					}
					if err := addFile(tw, filepath.Join(resultBase, "uniref"+suffix)); err != nil {
						return err
					}

					if taxonomy {
						if err := addFile(tw, filepath.Join(resultBase, "uniref_tax.tsv")); err != nil {
							return err
						}
					}

					if useTemplates {
						if err := addFile(tw, filepath.Join(resultBase, "pdb70.m8")); err != nil {
							return err
						}
					}

					if useEnv {
						if err := addFile(tw, filepath.Join(resultBase, "bfd.mgnify30.metaeuk30.smag30"+suffix)); err != nil {
							return err
						}
					}

					if err := addFile(tw, scriptPath); err != nil {
						return err
					}
				}
				*/
				
				return nil
			}()

			if err != nil {
				file.Close()
				return &JobExecutionError{err}
			}

			if err = file.Sync(); err != nil {
				file.Close()
				return &JobExecutionError{err}
			}

			if err = file.Close(); err != nil {
				return &JobExecutionError{err}
			}
		}

		if config.Verbose {
			log.Print("Process finished gracefully without error")
		}
		return nil
	case PairJob:
		resultBase := filepath.Join(config.Paths.Results, string(request.Id))

		scriptPath := filepath.Join(resultBase, "pair.sh")
		script, err := os.Create(scriptPath)
		if err != nil {
			return &JobExecutionError{err}
		}
		script.WriteString(`#!/bin/bash -e
MMSEQS="$1"
QUERY="$2"
BASE="$4"
DB1="$5"
CWD="$6"
SEARCH_PARAM="--num-iterations 3 --db-load-mode 2 -a --k-score 'seq:96,prof:80' -e 0.1 --max-seqs 10000"
EXPAND_PARAM="--expansion-mode 0 -e inf --expand-filter-clusters 0 --max-seq-id 0.95"
export MMSEQS_CALL_DEPTH=1
python3 "${CWD}/mmseqs-server/backend/get_intermediates.py" "${BASE}/job.fasta" /mnt/disks/colabfold-dbs/ColabFold/MsaServer/intermediate_store
"${MMSEQS}" createdb "${QUERY}" "${BASE}/qdb" --shuffle 0
"${MMSEQS}" pairaln "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign" "${BASE}/res_exp_realign_pair" --db-load-mode 2
"${MMSEQS}" align   "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign_pair" "${BASE}/res_exp_realign_pair_bt" --db-load-mode 2 -e inf -a
"${MMSEQS}" pairaln "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign_pair_bt" "${BASE}/res_final" --db-load-mode 2
"${MMSEQS}" convertalis "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign_pair_bt" "${BASE}/convertalis_tax" --format-output target,evalue,taxid,taxname,taxlineage --db-load-mode 2
"${MMSEQS}" convertalis "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_exp_realign_pair_bt" "${BASE}/convertalis_seq" --format-output target,tseq --db-load-mode 2
"${MMSEQS}" result2msa "${BASE}/qdb" "${DB1}.idx" "${BASE}/res_final" "${BASE}/pair.a3m" --db-load-mode 2 --msa-format-mode 6
"${MMSEQS}" unpackdb "${BASE}/pair.a3m" "${BASE}" --unpack-name-mode 0 --unpack-suffix .a3m
python3 "${CWD}/mmseqs-server/backend/convertalis_seq_to_tsv.py" "${BASE}/convertalis_seq" --pair --a3m_0 "${BASE}/0.a3m" --a3m_1 "${BASE}/1.a3m"
python3 mmseqs-server/backend/add_tax_to_msa.py "${BASE}/convertalis_tax" "${BASE}"
"${MMSEQS}" rmdb "${BASE}/qdb"
"${MMSEQS}" rmdb "${BASE}/qdb_h"
"${MMSEQS}" rmdb "${BASE}/res"
"${MMSEQS}" rmdb "${BASE}/res_exp"
"${MMSEQS}" rmdb "${BASE}/res_exp_realign"
"${MMSEQS}" rmdb "${BASE}/res_exp_realign_pair"
"${MMSEQS}" rmdb "${BASE}/res_exp_realign_pair_bt"
"${MMSEQS}" rmdb "${BASE}/res_final"
rm -rf -- "${BASE}/tmp"
`)
		err = script.Close()
		if err != nil {
			return &JobExecutionError{err}
		}
		mydir, err := os.Getwd()

		parameters := []string{
			"/bin/sh",
			scriptPath,
			config.Paths.Mmseqs,
			filepath.Join(resultBase, "job.fasta"),
			config.Paths.Databases,
			resultBase,
			config.Paths.ColabFold.Uniref,
			mydir,
		}

		cmd, done, err := execCommand(config.Verbose, parameters...)
		if err != nil {
			return &JobExecutionError{err}
		}

		select {
		case <-time.After(1 * time.Hour):
			if err := KillCommand(cmd); err != nil {
				log.Printf("Failed to kill: %s\n", err)
			}
			return &JobTimeoutError{}
		case err := <-done:
			if err != nil {
				return &JobExecutionError{err}
			}

			path := filepath.Join(filepath.Clean(config.Paths.Results), string(request.Id))
			file, err := os.Create(filepath.Join(path, "mmseqs_results_"+string(request.Id)+".tar.gz"))
			if err != nil {
				return &JobExecutionError{err}
			}

			err = func() (err error) {
				gw := gzip.NewWriter(file)
				defer func() {
					cerr := gw.Close()
					if err == nil {
						err = cerr
					}
				}()
				tw := tar.NewWriter(gw)
				defer func() {
					cerr := tw.Close()
					if err == nil {
						err = cerr
					}
				}()

				if err := addFile(tw, filepath.Join(resultBase, "0.a3m.tax")); err != nil {
					return err
				}

				if err := addFile(tw, filepath.Join(resultBase, "1.a3m.tax")); err != nil {
					return err
				}

				if err := addFile(tw, filepath.Join(resultBase, "convertalis_seq.tsv")); err != nil {
					return err
				}

				return nil
			}()

			if err != nil {
				file.Close()
				return &JobExecutionError{err}
			}

			if err = file.Sync(); err != nil {
				file.Close()
				return &JobExecutionError{err}
			}

			if err = file.Close(); err != nil {
				return &JobExecutionError{err}
			}
		}
		if config.Verbose {
			log.Print("Process finished gracefully without error")
		}
		return nil
	case IndexJob:
		file := filepath.Join(config.Paths.Databases, job.Path)
		params, err := ReadParams(file + ".params")
		if err != nil {
			return &JobExecutionError{err}
		}
		params.Status = StatusRunning
		err = SaveParams(file+".params", params)
		if err != nil {
			return &JobExecutionError{err}
		}
		err = CheckDatabase(file, params, config)
		if err != nil {
			params.Status = StatusError
			SaveParams(file+".params", params)
			return &JobExecutionError{err}
		}
		if config.Verbose {
			log.Println("Process finished gracefully without error")
		}
		params.Status = StatusComplete
		err = SaveParams(file+".params", params)
		if err != nil {
			return &JobExecutionError{err}
		}
		return nil
	default:
		return &JobInvalidError{}
	}
}

func worker(jobsystem JobSystem, config ConfigRoot) {
	log.Println("MMseqs2 worker")
	mailer := MailTransport(NullTransport{})
	if config.Mail.Mailer != nil {
		log.Println("Using " + config.Mail.Mailer.Type + " mail transport")
		mailer = config.Mail.Mailer.GetTransport()
	}

	var shouldExit int32 = 0
	if config.Worker.GracefulExit {
		go func() {
			sig := make(chan os.Signal, 1)
			signal.Notify(sig, syscall.SIGINT, syscall.SIGTERM)
			defer signal.Stop(sig)
			<-sig
			atomic.StoreInt32(&shouldExit, 1)
		}()
	}

	for {
		if config.Worker.GracefulExit && atomic.LoadInt32(&shouldExit) == 1 {
			return
		}
		ticket, err := jobsystem.Dequeue()
		if err != nil {
			if ticket != nil {
				log.Print(err)
			}
			time.Sleep(100 * time.Millisecond)
			continue
		}

		if ticket == nil && err == nil {
			time.Sleep(100 * time.Millisecond)
			continue
		}

		jobFile := filepath.Join(config.Paths.Results, string(ticket.Id), "job.json")

		f, err := os.Open(jobFile)
		if err != nil {
			jobsystem.SetStatus(ticket.Id, StatusError)
			log.Print(err)
			continue
		}

		var job JobRequest
		dec := json.NewDecoder(bufio.NewReader(f))
		err = dec.Decode(&job)
		f.Close()
		if err != nil {
			jobsystem.SetStatus(ticket.Id, StatusError)
			log.Print(err)
			continue
		}

		jobsystem.SetStatus(ticket.Id, StatusRunning)
		err = RunJob(job, config)
		mailTemplate := config.Mail.Templates.Success
		switch err.(type) {
		case *JobExecutionError, *JobInvalidError:
			jobsystem.SetStatus(ticket.Id, StatusError)
			log.Print(err)
			mailTemplate = config.Mail.Templates.Error
		case *JobTimeoutError:
			jobsystem.SetStatus(ticket.Id, StatusError)
			log.Print(err)
			mailTemplate = config.Mail.Templates.Timeout
		case nil:
			jobsystem.SetStatus(ticket.Id, StatusComplete)
		}
		if job.Email != "" {
			err = mailer.Send(Mail{
				config.Mail.Sender,
				job.Email,
				fmt.Sprintf(mailTemplate.Subject, string(ticket.Id)),
				fmt.Sprintf(mailTemplate.Body, string(ticket.Id)),
			})
			if err != nil {
				log.Print(err)
			}
		}
	}
}
