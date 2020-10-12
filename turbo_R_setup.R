options(java.parameters = "-Xmx6g")

# see also https://jangorecki.gitlab.io/data.cube/library/data.table/html/dcast.data.table.html
library(config)
library(dplyr)
library(ggplot2)
library(httr)
library(igraph)
library(jsonlite)
library(randomForest)
library(rdflib)
library(readr)
library(readxl)
library(reshape2)
library(RJDBC)
library(solrium)
library(ssh)
library(stringdist)
library(stringr)
library(tm)
library(uuid)

# train
library(splitstackshape)

### validation
library(ROCR)
library(caret)

# library(xgboost)
# # also try party and xgboot

# ensure that large integers aren't casted to scientific notation
#  for example when being posted into a SQL query
options(scipen = 999)

# make sure this is being read from the intended folder
# user's home?
# current working directory?

print("Default file path set to:")
print(getwd())

# what's wrong with these lines?
#Error in source("https://raw.githubusercontent.com/PennTURBO/turbo-globals/master/turbo_R_setup.R") : 
#  https://raw.githubusercontent.com/PennTURBO/turbo-globals/master/turbo_R_setup.R:44:18: unexpected input
#43: 
#44: pre_commit_tags <–

pre_commit_tags = readLines("../pre_commit_tags.txt")
pre_commit_status = readLines("../pre_commit_status.txt")

config.file <- "config/turbo_R_setup.yaml"

config <- config::get(file = config.file)

####

execution.timestamp <- as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S")
execution.timestamp <- strftime(execution.timestamp , "%Y-%m-%dT%H:%M:%SZ")

####

chunk.vec <- function(vec, chunk.count) {
  split(vec, cut(seq_along(vec), chunk.count, labels = FALSE))
}

make.table.frame <- function(my.vector) {
  temp <- table(my.vector)
  temp <- cbind.data.frame(names(temp), as.numeric(temp))
  colnames(temp) <- c('value', 'count')
  temp$value <- as.character(temp$value)
  return(temp)
}

label.table <- function() {
    temp <- table(term.label)
    temp <-
      cbind.data.frame(names(temp), as.numeric(temp))
    colnames(temp) <- c("label", "count")
    table(temp$count)
}

approximateTerm <- function(med.string) {
  params <- list(term = med.string, maxEntries = 50)
  r <-
    httr::GET(
      paste0("http://",
             rxnav.api.address,
             ":",
             rxnav.api.port,
             "/"),
      path = "REST/approximateTerm.json",
      query = params
    )
  r <- rawToChar(r$content)
  r <- jsonlite::fromJSON(r)
  r <- r$approximateGroup$candidate
  if (is.data.frame(r)) {
    r$query <- med.string
    Sys.sleep(0.1)
    return(r)
  }
}

bulk.approximateTerm <-
  function(strs = c("tylenol", "cisplatin", "benadryl", "rogaine")) {
    temp <- lapply(strs, function(current.query) {
      print(current.query)
      params <- list(term = current.query, maxEntries = 50)
      r <-
        httr::GET("http://localhost:4000/",
                  path = "REST/approximateTerm.json",
                  query = params)
      r <- rawToChar(r$content)
      r <- jsonlite::fromJSON(r)
      r <- r$approximateGroup$candidate
      if (is.data.frame(r)) {
        r$query <- current.query
        return(r)
      }
    })
    temp <-
      do.call(rbind.data.frame, temp)
    
    temp$rank <-
      as.numeric(as.character(temp$rank))
    temp$score <-
      as.numeric(as.character(temp$score))
    temp$rxcui <-
      as.numeric(as.character(temp$rxcui))
    temp$rxaui <-
      as.numeric(as.character(temp$rxaui))
    
    approximate.rxcui.tab <- table(temp$rxcui)
    approximate.rxcui.tab <-
      cbind.data.frame(names(approximate.rxcui.tab),
                       as.numeric(approximate.rxcui.tab))
    names(approximate.rxcui.tab) <- c("rxcui", "rxcui.count")
    approximate.rxcui.tab$rxcui <-
      as.numeric(as.character(approximate.rxcui.tab$rxcui))
    approximate.rxcui.tab$rxcui.freq <-
      approximate.rxcui.tab$rxcui.count / (sum(approximate.rxcui.tab$rxcui.count))
    
    
    approximate.rxaui.tab <- table(temp$rxaui)
    approximate.rxaui.tab <-
      cbind.data.frame(names(approximate.rxaui.tab),
                       as.numeric(approximate.rxaui.tab))
    names(approximate.rxaui.tab) <- c("rxaui", "rxaui.count")
    approximate.rxaui.tab$rxaui <-
      as.numeric(as.character(approximate.rxaui.tab$rxaui))
    approximate.rxaui.tab$rxaui.freq <-
      approximate.rxaui.tab$rxaui.count / (sum(approximate.rxaui.tab$rxaui.count))
    
    temp <-
      base::merge(x = temp, y = approximate.rxcui.tab)
    
    temp <-
      base::merge(x = temp, y = approximate.rxaui.tab)
    
    return(temp)
  }

bulk.rxaui.asserted.strings <-
  function(rxauis, chunk.count = rxaui.asserted.strings.chunk.count) {
    rxn.chunks <-
      chunk.vec(sort(unique(rxauis)), chunk.count)
    
    rxaui.asserted.strings <-
      lapply(names(rxn.chunks), function(current.index) {
        current.chunk <- rxn.chunks[[current.index]]
        tidied.chunk <-
          paste0("'", current.chunk, "'", collapse = ", ")
        
        rxnav.rxaui.strings.query <-
          paste0(
            "SELECT RXCUI as rxcui,
            RXAUI as rxaui,
            SAB ,
            SUPPRESS ,
            TTY ,
            STR
            from
            rxnorm_current.RXNCONSO r where RXAUI in ( ",
            tidied.chunk,
            ")"
          )
        
        temp <- dbGetQuery(rxnCon, rxnav.rxaui.strings.query)
        return(temp)
      })
    
    rxaui.asserted.strings <-
      do.call(rbind.data.frame, rxaui.asserted.strings)
    
    rxaui.asserted.strings[, c("rxcui", "rxaui")] <-
      lapply(rxaui.asserted.strings[, c("rxcui", "rxaui")],  as.numeric)
    
    rxaui.asserted.strings$STR.lc <-
      tolower(rxaui.asserted.strings$STR)
    
    return(rxaui.asserted.strings)
  }

get.string.dist.mat <- function(two.string.cols) {
  two.string.cols <- as.data.frame(two.string.cols)
  unique.string.combos <- unique(two.string.cols)
  distance.cols = c("lv", "lcs", "qgram", "cosine", "jaccard", "jw")
  distances <- lapply(distance.cols, function(one.meth) {
    print(one.meth)
    temp <-
      stringdist(
        a = two.string.cols[, 1],
        b = two.string.cols[, 2],
        method = one.meth,
        nthread = 4
      )
    return(temp)
  })
  distances <- do.call(cbind.data.frame, distances)
  colnames(distances) <- distance.cols
  two.string.cols <-
    cbind.data.frame(two.string.cols, distances)
  return(two.string.cols)
}

instantiate.and.upload <- function(current.task) {
  print(current.task)
  
  # more.specific <-
  #   config::get(file = "rxnav_med_mapping.yaml", config = current.task)
  
  more.specific <-
    config::get(file = config.file, config = current.task)
  
  predlist <- colnames(body[2:ncol(body)])
  print(predlist)
  
  current.model.rdf <- rdflib::rdf()
  
  placeholder <-
    apply(
      X = body,
      MARGIN = 1,
      FUN = function(current_row) {
        innerph <- lapply(predlist, function(current.pred) {
          rdflib::rdf_add(
            rdf = current.model.rdf,
            subject = current_row[[1]],
            predicate = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
            object = more.specific$my.class
          )
          temp <- current_row[[current.pred]]
          if (nchar(temp) > 0) {
            # print(paste0(current.pred, ':', temp))
            if (current.pred %in% more.specific$my.numericals) {
              temp <- as.numeric(temp)
            }
            rdflib::rdf_add(
              rdf = current.model.rdf,
              subject = current_row[[1]],
              predicate = paste0('http://example.com/resource/', current.pred),
              object = temp
            )
          }
        })
      }
    )
  
  rdf.file <- paste0(current.task, '.ttl')
  
  rdflib::rdf_serialize(rdf = current.model.rdf,
                        doc = rdf.file,
                        format = "turtle")
  
  post.dest <-
    paste0(
      more.specific$my.graphdb.base,
      '/repositories/',
      more.specific$my.selected.repo,
      '/rdf-graphs/service?graph=',
      URLencode(
        paste0('http://example.com/resource/',
               current.task),
        reserved = TRUE
      )
    )
  
  print(post.dest)
  
  post.resp <-
    httr::POST(
      url = post.dest,
      body = upload_file(rdf.file),
      content_type(more.specific$my.mappings.format),
      authenticate(
        more.specific$my.graphdb.username,
        more.specific$my.graphdb.pw,
        type = 'basic'
      )
    )
  
  print('Errors will be listed below:')
  print(rawToChar(post.resp$content))
  
}


import.from.local.file <-
  function(some.graph.name,
           some.local.file,
           some.rdf.format) {
    print(some.graph.name)
    print(some.local.file)
    print(some.rdf.format)
    post.dest <-
      paste0(
        config$my.graphdb.base,
        '/repositories/',
        config$my.selected.repo,
        '/rdf-graphs/service?graph=',
        some.graph.name
      )
    
    print(post.dest)
    
    post.resp <-
      httr::POST(
        url = post.dest,
        body = upload_file(some.local.file),
        content_type(some.rdf.format),
        authenticate(config$my.graphdb.username,
                     config$my.graphdb.pw,
                     type = 'basic')
      )
    
    print('Errors will be listed below:')
    print(rawToChar(post.resp$content))
  }

import.from.url <-   function(some.graph.name,
                              some.ontology.url,
                              some.rdf.format) {
  print(some.graph.name)
  print(some.ontology.url)
  print(some.rdf.format)
  
  if (nchar(some.rdf.format) > 0) {
    update.body <- paste0(
      '{
      "context": "',
      some.graph.name,
      '",
      "data": "',
      some.ontology.url,
      '",
      "format": "',
      some.rdf.format,
      '"
  }'
    )
  } else {
    update.body <- paste0('{
                          "context": "',
                          some.graph.name,
                          '",
                          "data": "',
                          some.ontology.url,
                          '"
  }')
  }
  
  cat("\n")
  cat(update.body)
  cat("\n\n")
  
  post.res <- POST(
    url.post.endpoint,
    body = update.body,
    content_type("application/json"),
    accept("application/json"),
    saved.authentication
  )
  
  cat(rawToChar(post.res$content))
  
}

get.context.report <- function() {
  context.report <- GET(
    url = paste0(
      config$my.graphdb.base,
      "/repositories/",
      config$my.selected.repo,
      "/contexts"
    ),
    saved.authentication
  )
  context.report <-
    jsonlite::fromJSON(rawToChar(context.report$content))
  context.report <-
    context.report$results$bindings$contextID$value
  return(context.report)
}

monitor.named.graphs <- function() {
  while (TRUE) {
    print(paste0(
      Sys.time(),
      ": '",
      last.post.status,
      "' submitted at ",
      last.post.time
    ))
    
    context.report <- get.context.report()
    
    pending.graphs <- sort(setdiff(expectation, context.report))
    
    # will this properly handle the case when the report is empty (NULL)?
    if (length(pending.graphs) == 0) {
      print("Update complete")
      break()
    }
    
    print(paste0("still waiting for: ", pending.graphs))
    
    print(paste0("Next check in ",
                 config$monitor.pause.seconds,
                 " seconds."))
    
    Sys.sleep(config$monitor.pause.seconds)
    
  }
}

q2j2df <-
  function(query,
           endpoint = config$my.graphdb.base,
           repo = config$my.selected.repo,
           auth = saved.authentication) {
    # query <- config$main.solr.query
    
    minquery <- gsub(pattern = " +",
                     replacement = " ",
                     x = query)
    
    rdfres <- httr::GET(
      url = paste0(endpoint,
                   "/repositories/",
                   repo),
      query = list(query = minquery),
      auth
    )
    
    # convert binary JSON SPARQL results to a minimal dataframe
    rdfres <-
      jsonlite::fromJSON(rawToChar(rdfres$content))
    rdfres <- rdfres$results$bindings
    
    
    rdfres <-
      do.call(what = cbind.data.frame, args = rdfres)
    keepers <- colnames(rdfres)
    keepers <- keepers[grepl(pattern = "value$", x = keepers)]
    rdfres <- rdfres[, keepers]
    
    # beautify column labels
    temp <-
      gsub(pattern = '\\.value$',
           replacement = '',
           x = colnames(rdfres))
    # temp <- gsub(pattern = '^.*\\$',
    #              replacement = '',
    #              x = temp)
    colnames(rdfres) <- temp
    
    return(rdfres)
    
  }

url.post.endpoint <-
  paste0(
    config$my.graphdb.base,
    "/rest/data/import/upload/",
    config$my.selected.repo,
    "/url"
  )

update.endpoint <-
  paste0(config$my.graphdb.base,
         "/repositories/",
         config$my.selected.repo,
         "/statements")

select.endpoint <-
  paste0(config$my.graphdb.base,
         "/repositories/",
         config$my.selected.repo)


####


saved.authentication <-
  authenticate(config$my.graphdb.username,
               config$my.graphdb.pw,
               type = "basic")

####


rxnDriver <-
  JDBC(driverClass = "com.mysql.cj.jdbc.Driver",
       classPath = config$mysql.jdbc.path)

# # i keep re-doing this thorugh other scripts
# rxnCon <-
#   dbConnect(
#     rxnDriver,
#     paste0(
#       "jdbc:mysql://",
#       config$rxnav.mysql.address,
#       ":",
#       config$rxnav.mysql.port
#     ),
#     config$rxnav.mysql.user,
#     config$rxnav.mysql.pw
#   )

####

### get mappings with BioPortal
# or string-search somewhere?
# start with public endpoint but eventually switch to appliance

#### these are functioning like globals so they don't have to be passed to the function
api.base.uri <- "http://data.bioontology.org/ontologies"
api.ontology.name <- "LOINC"
term.ontology.name <- "LNC"
term.base.uri <-
  paste0("http://purl.bioontology.org/ontology",
         "/",
         term.ontology.name)
api.family <- "classes"
# source.term <- "http://purl.bioontology.org/ontology/LNC/LP17698-9"
api.method <- "mappings"
# what are the chances that a mapping query will return 0 mappings, or that it will return multiple pages?

bp.map.retreive.and.parse <- function(term.list) {
  outer <- lapply(term.list, function(current.term) {
    # current.term <- "LP102314-4"
    # current.term <-"LP40488-6"
    # current.term <-"LP417915-8"
    
    print(current.term)
    current.uri <- paste0(term.base.uri, "/", current.term)
    encoded.term <- URLencode(current.uri, reserved = TRUE)
    prepared.get <-
      paste(api.base.uri,
            api.ontology.name,
            api.family,
            encoded.term,
            api.method,
            sep = "/")
    mapping.res.list <-
      httr::GET(url = prepared.get,
                add_headers(
                  Authorization = paste0("apikey token=", config$public.bioportal.api.key)
                ))
    
    print(mapping.res.list$status_code)
    
    if (mapping.res.list$status_code == 200) {
      mapping.res.list <- rawToChar(mapping.res.list$content)
      
      mapping.res.list <- jsonlite::fromJSON(mapping.res.list)
      
      # print(head(mapping.res.list))
      
      if (length(mapping.res.list) > 0) {
        # CUI, LOOM, "same URI", etc. Probably only LOOM will be useful
        mapping.methods <- mapping.res.list$source
        
        source.target.details <-
          lapply(mapping.res.list$classes, function(current.mapping) {
            source.target.terms <- current.mapping$`@id`
            source.target.ontologies <-
              current.mapping$links$ontology
            return(c(
              rbind(source.target.terms, source.target.ontologies)
            ))
          })
        
        source.target.details <-
          do.call(rbind.data.frame, source.target.details)
        colnames(source.target.details) <-
          c("source.term",
            "source.ontology",
            "target.term",
            "target.ontology")
        
        source.target.details <-
          cbind.data.frame(source.target.details, mapping.methods)
        return(source.target.details)
      }
    }
  })
}

bioportal.string.search <- function(current.string) {
  # current.string <- 'asthma'
  print(current.string)
  prepared.get <-
    paste0(
      'http://data.bioontology.org/search?q=',
      current.string  ,
      '&include=prefLabel,synonym',
      '&pagesize=999'
    )
  prepared.get <- URLencode(prepared.get, reserved = FALSE)
  search.res.list <-
    httr::GET(url = prepared.get,
              add_headers(
                Authorization = paste0("apikey token=", config$public.bioportal.api.key)
              ))
  
  search.res.list <- rawToChar(search.res.list$content)
  search.res.list <- jsonlite::fromJSON(search.res.list)
  search.res.list <- search.res.list$collection
  
  # print(search.res.list$links$ontology)
  
  if (is.data.frame(search.res.list)) {
    if (nrow(search.res.list) > 0) {
      ontology <- search.res.list$links$ontology
      #  , 'ontologyType'
      search.res.list <- search.res.list[, c('prefLabel', '@id')]
      colnames(search.res.list) <- c('prefLabel', 'iri')
      search.res.list <-
        cbind.data.frame(search.res.list, 'ontology' = ontology)
      search.res.list$rank <- 1:nrow(search.res.list)
      return(search.res.list)
    }
  }
}


# see https://www.ebi.ac.uk/ols/docs/api
ols.serch.term.labels.universal <-
  function(current.string,
           current.id,
           strip.final.s = FALSE,
           ontology.filter,
           kept.row.count = 9,
           req.exact = 'false') {
    if (strip.final.s) {
      current.string <-
        sub(pattern = "s$",
            replacement = "",
            x = current.string)
    }
    
    # singular.lc <- current.string
    
    # or just try url encoding?
    
    # substitute 'spp$' or 'sp$' with ''  for genus-level NCBI taxon entities
    # that porabialy isn't desirable in general
    # and should be really clear to users fo thsi function
    
    current.string <-
      gsub(pattern = " sp$",
           replacement = "",
           x = current.string)
    
    current.string <-
      gsub(pattern = " spp$",
           replacement = "",
           x = current.string)
    
    singular.lc <- current.string
    print(singular.lc)
    
    current.string <-
      gsub(pattern = "[[:punct:] ]",
           replacement = ",",
           x = current.string)
    
    
    print(current.string)
    
    prepared.query <- paste0(
      "https://www.ebi.ac.uk/ols/api/search?q={",
      current.string,
      "}&type=class&local=true",
      ontology.filter ,
      "&rows=",
      kept.row.count,
      '&exact=',
      req.exact,
      "&fieldList=iri,short_form,obo_id,ontology_name,ontology_prefix,label,synonym,annotations,annotations_trimmed",
      "&query_fields=label,synonym,annotations,annotations_trimmed"
    )
    
    # print(prepared.query)
    
    ols.attempt <-
      httr::GET(prepared.query)
    
    ols.attempt <- ols.attempt$content
    ols.attempt <- rawToChar(ols.attempt)
    ols.attempt <- jsonlite::fromJSON(ols.attempt)
    ols.attempt <- ols.attempt$response$docs
    if (is.data.frame(ols.attempt)) {
      if (nrow(ols.attempt) > 0) {
        ols.attempt$query <- singular.lc
        ols.attempt$loinc.part <- current.id
        
        ols.attempt$rank <- 1:nrow(ols.attempt)
        ols.attempt$label <- tolower(ols.attempt$label)
        ols.attempt$query <- tolower(ols.attempt$query)
        
        return(ols.attempt)
      }
    }
  }

# updates current.component.mapping.frame
update.accounting <- function(data,
                              loinc.part.code,
                              loinc.part.name ,
                              obo.uri,
                              obo.label,
                              rank,
                              justification) {
  print("before update")
  print(length(current.component.mapping.complete))
  print(length(current.needs.component.mapping))
  print(nrow(current.component.mapping.frame))
  
  print("update row count")
  print(nrow(data))
  
  bare.lpc <- unlist(data[, loinc.part.code])
  
  current.component.mapping.complete <<-
    union(current.component.mapping.complete, bare.lpc)
  
  current.needs.component.mapping <<-
    setdiff(current.needs.component.mapping, bare.lpc)
  
  matches.external.cols <-
    data[, c(loinc.part.code, loinc.part.name, obo.uri, obo.label, rank)]
  matches.external.cols$justification <- justification
  
  colnames(matches.external.cols) <-
    c(
      'loinc.part.code',
      'loinc.part.name',
      'obo.uri',
      'obo.label',
      'rank',
      'justification'
    )
  
  current.component.mapping.frame <<-
    rbind.data.frame(current.component.mapping.frame, matches.external.cols)
  
  print("after update")
  print(length(current.component.mapping.complete))
  print(length(current.needs.component.mapping))
  print(nrow(current.component.mapping.frame))
  
  print(sort(table(
    current.component.mapping.frame$justification
  )))
}

split.details <- function(PartTypeNameVal, acceptable.details) {
  has.details <-
    LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == PartTypeNameVal]
  
  has.details <-
    intersect(has.details, ehr.with.loinc.parts$LOINC)
  
  print(sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% has.details &
                                            LoincPartLink$PartTypeName == PartTypeNameVal])))
  
  acceptable.details.codes <-
    LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == PartTypeNameVal &
                                LoincPartLink$PartName %in% acceptable.details]
  
  unacceptable.details <-
    setdiff(has.details, acceptable.details.codes)
  
  ehr.with.loinc.parts <-
    ehr.with.loinc.parts[!ehr.with.loinc.parts$LOINC %in% unacceptable.details , ]
  
  print(nrow(ehr.with.loinc.parts))
  
  detail.frame <-
    unique(LoincPartLink[LoincPartLink$PartTypeName == PartTypeNameVal  &
                           LoincPartLink$LoincNumber %in%  ehr.with.loinc.parts$LOINC , c('LoincNumber', "PartName")])
  
  detail.prep <- detail.frame
  
  detail.prep$placeholder <- TRUE
  
  detail.cast <-
    dcast(data = detail.prep,
          formula = LoincNumber ~ PartName,
          value.var = 'placeholder')
  detail.cast$detail.count <-
    rowSums(detail.cast[,-1], na.rm = TRUE)
  
  # always.keep <- c('LoincNumber', 'detail.count')
  
  detail.followup.cols <-
    c(setdiff(colnames(detail.cast), acceptable.details))
  
  detail.followup <- detail.cast[, detail.followup.cols]
  
  detail.cast <-
    detail.cast[, union('LoincNumber', acceptable.details)]
  
  return(list(detail.cast = detail.cast, detail.followup = detail.followup))
  
}

# # fixme
# selected.columns <- divisors
# part.name <- 'COMPONENT'

tl.augmenter <- function(selected.columns, part.name) {
  details.frame <- pre.ready.for.robot[, selected.columns]
  
  details.key <- pre.ready.for.robot$LOINC
  details.tally <- rowSums(details.frame, na.rm = TRUE)
  details.tally <- cbind.data.frame(details.key, details.tally)
  
  details.frame <-
    cbind.data.frame(details.key, details.frame)
  
  details.melt <-
    melt(data = details.frame, id.vars = 'details.key')
  details.melt[] <- lapply(X = details.melt[], FUN = as.character)
  details.melt <-
    as.data.frame(details.melt[complete.cases(details.melt),])
  # table(details.melt$variable)
  # print(length(unique(details.melt$details.key)))
  dm.check <- make.table.frame(details.melt$details.key)
  dm.check <- dm.check$value[dm.check$count > 1]
  dm.singles <-
    details.melt[(!(details.melt$details.key %in% dm.check)) ,]
  dm.check <-
    details.melt[details.melt$details.key %in% dm.check ,]
  
  if (nrow(dm.check) > 0) {
    dm.check$nchar <- nchar(dm.check$variable)
    
    dm.longest <- aggregate(dm.check$nchar,
                            by = list(dm.check$details.key),
                            FUN = max)
    colnames(dm.longest) <- c('details.key', 'nchar')
    
    dm.check <-
      base::merge(x = dm.check , y = dm.longest)
    
    dm.singles <- dm.singles[, colnames(details.melt)]
    dm.check <- dm.check[, colnames(details.melt)]
    details.melt <- rbind.data.frame(dm.singles, dm.check)
  }
  
  rfr.min <-
    as.data.frame(pre.ready.for.robot[, c('LOINC', part.name)])
  # print(length(unique(rfr.min$LOINC)))
  # temp <- make.table.frame(rfr.min$LOINC)
  
  details.join <-
    left_join(x = rfr.min,
              y = details.melt,
              by = c("LOINC" = "details.key"))
  
  details.join <- details.join[order(details.join$LOINC), ]
  
  return(details.join$variable)
}

rxnCon <- NULL

# todo paramterize connection and query string
# how to user conenction parpatmeron LHS or assignment?
rxnav.test.and.refresh <- function() {
  local.q <- "select RSAB from rxnorm_current.RXNSAB r"
  tryCatch({
    dbGetQuery(rxnCon, local.q)
  }, warning = function(w) {
    
  }, error = function(e) {
    print(e)
    print("trying to reconnect")
    rxnCon <<- dbConnect(
      rxnDriver,
      paste0(
        "jdbc:mysql://",
        config$rxnav.mysql.address,
        ":",
        config$rxnav.mysql.port
      ),
      config$rxnav.mysql.user,
      config$rxnav.mysql.pw
    )
    dbGetQuery(rxnCon, local.q)
  }, finally = {
    
  })
}

build.source.med.classifications.annotations <-
  function(version.list,
           onto.iri,
           onto.file,
           onto.file.format) {

    # cat(config$source.med.classifications.onto.comment)
    
    annotation.model <- rdf()
    
    rdflib::rdf_add(
      rdf = annotation.model,
      subject = onto.iri,
      predicate = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
      object = "http://www.w3.org/2002/07/owl#Ontology"
    )
    
    rdflib::rdf_add(
      rdf = annotation.model,
      subject = onto.iri,
      predicate = "http://purl.org/dc/terms/created",
      object = version.list$created
    )
    
    rdflib::rdf_add(
      rdf = annotation.model,
      subject = onto.iri,
      predicate = "http://www.w3.org/2002/07/owl#versionInfo",
      object = version.list$versioninfo
    )
    
    rdflib::rdf_add(
      rdf = annotation.model,
      subject = onto.iri,
      predicate = "http://www.w3.org/2000/01/rdf-schema#comment",
      object = config$source.med.classifications.onto.comment
    )
    
    rdf_serialize(rdf = annotation.model,
                  doc = onto.file,
                  format = onto.file.format)
  }

